#' Derive Death Cause
#'
#' Derive death cause (`DTHCAUS`) and add traceability variables if required.
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#' @param ... Objects of class "dthcaus_source" created by [`dthcaus_source()`].
#'
#' @details
#' This function derives `DTHCAUS` along with the user-defined traceability
#' variables, if required. If a subject has death info from multiple sources,
#' the one from the source with the earliest death date will be used. If dates are
#' equivalent, the first source will be kept, so the user should provide the inputs in
#' the preferred order.
#'
#' @keywords adsl
#'
#' @author
#' Shimeng Huang
#' Samia Kabi
#'
#' @return The input dataset with `DTHCAUS` variable added.
#'
#' @export
#'
#' @seealso [dthcaus_source()]
#'
#' @examples
#' adsl <- tibble::tribble(
#'   ~STUDYID, ~USUBJID,
#'   "STUDY01", "PAT01",
#'   "STUDY01", "PAT02"
#' )
#' ae <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
#'   "STUDY01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
#' )
#' ds <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01"
#' )
#'
#' # Derive `DTHCAUS` only
#' src_ae <- dthcaus_source(
#'   dataset = ae,
#'   filter = AEOUT == "FATAL",
#'   date_var = AEDTHDTC,
#'   mode = "first",
#'   dthcaus = AEDECOD
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset = ds,
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDTC,
#'   mode = "first",
#'   dthcaus = DSTERM
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds)
#'
#' # Derive `DTHCAUS` and add traceability variables
#' src_ae <- dthcaus_source(
#'   dataset = ae,
#'   filter = AEOUT == "FATAL",
#'   date_var = AEDTHDTC,
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset = ds,
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date_var = DSSTDTC,
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   traceabilty_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds)
derive_var_dthcaus <- function(dataset, ...) {
  sources <- list(...)
  walk(sources, validate_dthcaus_source)

  warn_if_vars_exist(dataset, "DTHCAUS")

  # process each source
  add_data <- vector("list", length(sources))
  for (ii in seq_along(sources)) {
    if (!is.null(sources[[ii]]$filter)) {
      add_data[[ii]] <- sources[[ii]]$dataset %>%
        filter(!!sources[[ii]]$filter)
    } else {
      add_data[[ii]] <- sources[[ii]]$dataset
    }

    # if several death records, use the first/last according to 'mode'
    add_data[[ii]] <- add_data[[ii]] %>%
      filter_extreme(
        order = vars(!!sources[[ii]]$date),
        by_vars = vars(USUBJID),
        mode = sources[[ii]]$mode
      )

    add_data[[ii]] <- add_data[[ii]] %>%
      mutate(
        temp_source_nr = ii,
        temp_date = !!sources[[ii]]$date,
        DTHCAUS = !!sources[[ii]]$dthcaus
      )

    # add traceabilty param if required
    # inconsitent traceability lists issue a warning
    if (ii > 1) {
      warn_if_inconsistent_list(
        base = sources[[ii - 1]]$traceabilty,
        compare = sources[[ii]]$traceabilty,
        list_name = "dthcaus_source()",
        i = ii
      )
    }
    if (!is.null(sources[[ii]]$traceabilty)) {
      warn_if_vars_exist(dataset, names(sources[[ii]]$traceabilty))
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(
          USUBJID,
          temp_source_nr,
          temp_date,
          DTHCAUS,
          !!!sources[[ii]]$traceabilty
        )
    }
    else {
      add_data[[ii]] <- add_data[[ii]] %>%
        select(USUBJID, temp_source_nr, temp_date, DTHCAUS)
    }
  }
  # if a subject has multiple death info, keep the one from the first source
  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(
      order = vars(temp_date, temp_source_nr),
      by_vars = vars(USUBJID),
      mode = "first"
    ) %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}

#' Create an `dthcaus_source` object
#'
#' @param dataset A data.frame containing a source dataset.
#' @param filter An expression used for filtering `dataset`.
#' @param date_var A character vector to be used for sorting `dataset`.
#' @param mode One of `"first"` or `"last"`.
#' @param dthcaus A variable name or a string literal --- if a variable name, e.g., `AEDECOD`,
#'   it is the variable in the source dataset to be used to assign values to
#'   `DTHCAUS`; if a string literal, e.g. `"Adverse Event"`, it is the fixed value
#'   to be assigned to `DTHCAUS`.
#' @param traceabilty_vars A named list returned by [`vars()`] listing the traceability
#'  variables, e.g. `vars(DTHDOM = "DS", DTHSEQ = DSSEQ)`.
#'
#' @author Shimeng Huang
#'
#' @seealso [`derive_var_dthcaus()`]
#'
#' @export
#'
#' @return An object of class "dthcaus_source".
dthcaus_source <- function(dataset,
                           filter,
                           date_var,
                           mode = "first",
                           dthcaus,
                           traceabilty_vars = NULL) {
  out <- list(
    dataset = dataset,
    filter = enquo(filter),
    date = enquo(date_var),
    mode = mode,
    dthcaus = enquo(dthcaus),
    traceabilty = traceabilty_vars
  )
  class(out) <- c("dthcaus_source", "list")
  validate_dthcaus_source(out)
}

#' Validate an object is indeed a `dthcaus_source` object
#'
#' @param x An object to be validated.
#'
#' @author Shimeng Huang
#'
#' @return The original object.
validate_dthcaus_source <- function(x) {
  assert_that(inherits(x, "dthcaus_source"))
  values <- unclass(x)
  assert_that(is.data.frame(values$dataset))
  assert_that(is_expr(values$filter))
  assert_that(is_expr(values$date))
  assert_that(values$mode %in% c("first", "last"))
  assert_that(quo_is_symbol(values$dthcaus) | is.character(quo_get_expr(values$dthcaus)))
  assert_that(is.list(values$traceabilty) | is.null(values$traceability))
  x
}
