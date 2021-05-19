#' Derive Death Cause and Source Domain
#'
#' Derives death cause (`DTHCAUS`) and source domain (`DTHDOM`).
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#'
#' @param ... Objects of class "dthcaus_source" created by `dthcaus_source`.
#'
#' Note: If a subject has death info from multiple sources, the one from the
#'   first source will be kept.
#'
#' @keywords adsl
#'
#' @author Shimeng Huang
#'
#' @return The input dataset with `DTHCAUS` and `DTHDOM` added.
#'
#' @export
#'
#' @seealso [dthcaus_source()]
#'
#' @examples
#' adsl <- tibble::tribble(
#'  ~STUDYID, ~USUBJID,
#'  "STUDY01", "PAT01",
#'  "STUDY01", "PAT02"
#' )
#' ae <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
#'   "STUDY01", "PAT01", "SUDDEN DEATH" , "FATAL", "2021-04-04"
#' )
#' ds <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01"
#' )
#'
#' src_ae <- dthcaus_source(
#'   dataset = ae,
#'   filter = expr(AEOUT == "FATAL"),
#'   order = exprs(AEDTHDTC),
#'   mode = "first",
#'   dthdom = "AE",
#'   dthcaus = expr(AEDECOD)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset = ds,
#'   filter = expr(DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM)),
#'   order = exprs(DSSTDTC),
#'   mode = "first",
#'   dthdom = "DS",
#'   dthcaus = expr(DSTERM)
#' )
#'
#' derive_dthcaus(adsl, src_ae, src_ds)
derive_dthcaus <- function(dataset, ...) {

  sources <- list(...)

  # check and clean up input
  for (ss in sources) {
    validate_dthcaus_source(ss)
  }

  # process each source
  add_data <- vector("list", length(sources))
  for (ii in seq_along(sources)) {

    if (!is.null(sources[[ii]]$filter)) {
      add_data[[ii]] <- sources[[ii]]$dataset %>%
        filter(!!!(sources[[ii]]$filter))
    }
    else {
      add_data[[ii]] <- sources[[ii]]$dataset
    }

    if (!is.null(sources[[ii]]$order)){
      add_data[[ii]] <- add_data[[ii]] %>%
        filter_extreme(order = sources[[ii]]$order,
                       by_vars = exprs(USUBJID),
                       mode = sources[[ii]]$mode)
    }

    if (typeof(sources[[ii]]$dthcaus) == "symbol") {
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(USUBJID,
                  temp_source_nr = ii,
                  DTHDOM = sources[[ii]]$dthdom,
                  DTHCAUS := !!sources[[ii]]$dthcaus)
    }
    else {
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(USUBJID,
                  temp_source_nr = ii,
                  DTHDOM = sources[[ii]]$dthdom,
                  DTHCAUS = sources[[ii]]$dthcaus)
    }
  }

  # if a subject has multiple death info, keep the one from the first source
  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(order = exprs(temp_source_nr),
                   by_vars = exprs(USUBJID),
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}

#' Create an `dthcaus_source` object
#'
#' @param dataset A data.frame containint a source dataset.
#' @param filter A symbol returned by `expr` to be used for filtering `dataset`.
#' @param order Alist returned by `exprs` to be used for sorting `dataset`.
#' @param mode One of "first" or "last".
#' @param dthdom Name of the source domain as a string, e.g. "AE".
#' @param dthcaus A symbol or a string --- if a symbol, e.g., `expr(AEDECOD)`,
#'   it is the variable in the source dataset to be used to assign values to
#'   `DTHCAUS`; if a string, it is the fixed value to be assigned to `DTHCAUS`.
#'
#' @author Shimeng Huang
#'
#' @return An object of class "dthcaus_source".
dthcaus_source <- function(dataset, filter, order, mode = "first",
                           dthdom, dthcaus) {
  out <- list(
    dataset = dataset,
    filter = filter,
    order = order,
    mode = mode,
    dthdom = dthdom,
    dthcaus = dthcaus
  )
  class(out) <- "dthcaus_source"
  validate_dthcaus_source(out)
}

validate_dthcaus_source <- function(x) {
  stopifnot(inherits(x, "dthcaus_source"))
  values <- unclass(x)
  stopifnot("data.frame" %in% class(values$dataset))
  stopifnot(is_expr(values$filter))
  stopifnot(is_unnamed_exprs(values$order))
  match.arg(values$mode, c("first", "last"))
  stopifnot(is.character(values$dthdom))
  stopifnot(!is.symbol(values$dthcaus) | !is.character(values$dthcaus))
  x
}
