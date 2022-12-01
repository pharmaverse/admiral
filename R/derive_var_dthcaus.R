#' Derive Death Cause
#'
#' Derive death cause (`DTHCAUS`) and add traceability variables if required.
#'
#' @param dataset Input dataset.
#'
#'   The variables specified by `subject_keys` are required.
#'
#' @param source_datasets A named `list` containing datasets in which to search for the
#'   death cause
#'
#' @param ... Objects of class "dthcaus_source" created by [`dthcaus_source()`].
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @details
#' This function derives `DTHCAUS` along with the user-defined traceability
#' variables, if required. If a subject has death info from multiple sources,
#' the one from the source with the earliest death date will be used. If dates are
#' equivalent, the first source will be kept, so the user should provide the inputs in
#' the preferred order.
#'
#' @family der_adsl
#' @keywords der_adsl
#'
#' @author
#' Shimeng Huang, Samia Kabi, Thomas Neitmann, Tamara Senior
#'
#' @return The input dataset with `DTHCAUS` variable added.
#'
#' @export
#'
#' @seealso [dthcaus_source()]
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adsl <- tribble(
#'   ~STUDYID,  ~USUBJID,
#'   "STUDY01", "PAT01",
#'   "STUDY01", "PAT02",
#'   "STUDY01", "PAT03"
#' )
#' ae <- tribble(
#'   ~STUDYID,  ~USUBJID, ~AESEQ, ~AEDECOD,       ~AEOUT,  ~AEDTHDTC,
#'   "STUDY01", "PAT01",  12,     "SUDDEN DEATH", "FATAL", "2021-04-04"
#' ) %>%
#'   mutate(
#'     AEDTHDT = ymd(AEDTHDTC)
#'   )
#' ds <- tribble(
#'   ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
#'   "STUDY01", "PAT03", 1, "DEATH", "POST STUDY REPORTING OF DEATH", "2022-03-03"
#' ) %>%
#'   mutate(
#'     DSSTDT = ymd(DSSTDTC)
#'   )
#'
#' # Derive `DTHCAUS` only - for on-study deaths only
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = AEDTHDT,
#'   mode = "first",
#'   dthcaus = AEDECOD
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = DSTERM
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))
#'
#' # Derive `DTHCAUS` and add traceability variables - for on-study deaths only
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = AEDTHDT,
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))
#'
#' # Derive `DTHCAUS` as above - now including post-study deaths with different `DTHCAUS` value
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = AEDTHDT,
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' src_ds_post <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & DSTERM == "POST STUDY REPORTING OF DEATH",
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = "POST STUDY: UNKNOWN CAUSE",
#'   traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds, src_ds_post, source_datasets = list(ae = ae, ds = ds))
derive_var_dthcaus <- function(dataset,
                               ...,
                               source_datasets,
                               subject_keys = get_admiral_option("subject_keys")) {
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)
  assert_list_of(source_datasets, "data.frame")
  sources <- list(...)
  assert_list_of(sources, "dthcaus_source")

  source_names <- names(source_datasets)
  assert_list_element(
    list = sources,
    element = "dataset_name",
    condition = dataset_name %in% source_names,
    source_names = source_names,
    message_text = paste0(
      "The dataset names must be included in the list specified for the ",
      "`source_datasets` parameter.\n",
      "Following names were provided by `source_datasets`:\n",
      enumerate(source_names, quote_fun = squote)
    )
  )

  warn_if_vars_exist(dataset, "DTHCAUS")

  # process each source
  add_data <- vector("list", length(sources))
  for (ii in seq_along(sources)) {
    source_dataset_name <- sources[[ii]]$dataset_name
    source_dataset <- source_datasets[[source_dataset_name]]
    add_data[[ii]] <- source_dataset %>%
      filter_if(sources[[ii]]$filter)

    assert_date_var(
      dataset = add_data[[ii]],
      var = !!sources[[ii]]$date,
      dataset_name = source_dataset_name
    )

    # if several death records, use the first/last according to 'mode'
    tmp_source_nr <- get_new_tmp_var(dataset)
    tmp_date <- get_new_tmp_var(dataset)
    add_data[[ii]] <- add_data[[ii]] %>%
      filter_extreme(
        order = vars(!!sources[[ii]]$date, !!!sources[[ii]]$order),
        by_vars = subject_keys,
        mode = sources[[ii]]$mode
      ) %>%
      mutate(
        !!tmp_source_nr := ii,
        !!tmp_date := !!sources[[ii]]$date,
        DTHCAUS = !!sources[[ii]]$dthcaus
      )

    # add traceability param if required
    # inconsistent traceability lists issue a warning
    if (ii > 1) {
      warn_if_inconsistent_list(
        base = sources[[ii - 1]]$traceability,
        compare = sources[[ii]]$traceability,
        list_name = "dthcaus_source()",
        i = ii
      )
    }
    if (!is.null(sources[[ii]]$traceability)) {
      warn_if_vars_exist(source_dataset, names(sources[[ii]]$traceability))
      assert_data_frame(source_dataset, required_vars = get_source_vars(sources[[ii]]$traceability))
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(
          !!!subject_keys,
          !!tmp_source_nr,
          !!tmp_date,
          DTHCAUS,
          !!!sources[[ii]]$traceability
        )
    } else {
      add_data[[ii]] <- add_data[[ii]] %>%
        select(!!!subject_keys, !!tmp_source_nr, !!tmp_date, DTHCAUS)
    }
  }
  # if a subject has multiple death info, keep the one from the first source
  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(
      order = vars(!!tmp_date, !!tmp_source_nr),
      by_vars = subject_keys,
      mode = "first"
    ) %>%
    remove_tmp_vars()

  derive_vars_merged(dataset, dataset_add = dataset_add, by_vars = subject_keys)
}

#' Create a `dthcaus_source` Object
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the death cause.
#'
#' @param filter An expression used for filtering `dataset`.
#'
#' @param date A date or datetime variable to be used for sorting `dataset`.
#'
#' @param order Sort order
#'
#'   Additional variables to be used for sorting the `dataset` which is ordered by the
#'   `date` and `order`. Can be used to avoid duplicate record warning.
#'
#'   *Default*: `NULL`
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'   created by `vars()`, e.g., `vars(ADT, desc(AVAL))` or `NULL`
#'
#' @param mode One of `"first"` or `"last"`.
#' Either the `"first"` or `"last"` observation is preserved from the `dataset`
#' which is ordered by `date`.
#'
#' @param dthcaus A variable name or a string literal --- if a variable name, e.g., `AEDECOD`,
#'   it is the variable in the source dataset to be used to assign values to
#'   `DTHCAUS`; if a string literal, e.g. `"Adverse Event"`, it is the fixed value
#'   to be assigned to `DTHCAUS`.
#'
#' @param traceability_vars A named list returned by [`vars()`] listing the traceability variables,
#' e.g. `vars(DTHDOM = "DS", DTHSEQ = DSSEQ)`.
#' The left-hand side (names of the list elements) gives the names of the traceability variables
#' in the returned dataset.
#' The right-hand side (values of the list elements) gives the values of the traceability variables
#' in the returned dataset.
#' These can be either strings or symbols referring to existing variables.
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @author Shimeng Huang
#'
#' @export
#'
#' @seealso [derive_var_dthcaus()]
#'
#' @return An object of class "dthcaus_source".
#'
#' @examples
#' # Deaths sourced from AE
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = AEDTHDT,
#'   mode = "first",
#'   dthcaus = AEDECOD
#' )
#'
#' # Deaths sourced from DS
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH",
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = DSTERM
#' )
dthcaus_source <- function(dataset_name,
                           filter,
                           date,
                           order = NULL,
                           mode = "first",
                           dthcaus,
                           traceability_vars = NULL) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    order = assert_order_vars(order, optional = TRUE),
    mode = assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE),
    dthcaus = assert_symbol(enquo(dthcaus)) %or% assert_character_scalar(dthcaus),
    traceability = assert_varval_list(traceability_vars, optional = TRUE)
  )
  class(out) <- c("dthcaus_source", "source", "list")
  out
}
