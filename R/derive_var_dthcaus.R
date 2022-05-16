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
#' @keywords derivation adsl
#'
#' @author
#' Shimeng Huang, Samia Kabi, Thomas Neitmann
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
#'   "STUDY01", "PAT02",
#'   "STUDY01", "PAT03"
#' )
#' ae <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
#'   "STUDY01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
#' )
#' ds <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
#'   "STUDY01", "PAT03", 1, "DEATH", "POST STUDY REPORTING OF DEATH", "2022-03-03"
#' )
#'
#' # Derive `DTHCAUS` only - for on-study deaths only
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = AEDTHDTC,
#'   mode = "first",
#'   dthcaus = AEDECOD
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDTC,
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
#'   date = AEDTHDTC,
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDTC,
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
#'   date = AEDTHDTC,
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDTC,
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' src_ds_post <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & DSTERM == "POST STUDY REPORTING OF DEATH",
#'   date = DSSTDTC,
#'   mode = "first",
#'   dthcaus = "POST STUDY: UNKNOWN CAUSE",
#'   traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds, src_ds_post, source_datasets = list(ae = ae, ds = ds))
derive_var_dthcaus <- function(dataset,
                               ...,
                               source_datasets,
                               subject_keys = vars(STUDYID, USUBJID)) {
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
    if (!is.null(sources[[ii]]$date_imputation) & !quo_is_null(sources[[ii]]$filter)) {
      add_data[[ii]] <- source_dataset %>%
        filter(!!sources[[ii]]$filter) %>%
        mutate(temp_date = impute_dtc(dtc = !!sources[[ii]]$date,
                                      date_imputation = !!sources[[ii]]$date_imputation))
    } else if (!is.null(sources[[ii]]$date_imputation)) {
      add_data[[ii]] <- source_dataset %>%
        mutate(temp_date = impute_dtc(dtc = !!sources[[ii]]$date,
                                      date_imputation = !!sources[[ii]]$date_imputation))
    } else if (!quo_is_null(sources[[ii]]$filter)) {
      add_data[[ii]] <- source_dataset %>%
        filter(!!sources[[ii]]$filter) %>%
        mutate(temp_date = !!sources[[ii]]$date)
    } else {
      add_data[[ii]] <- source_dataset %>%
        mutate(temp_date = !!sources[[ii]]$date)
    }

    # if several death records, use the first/last according to 'mode'
    add_data[[ii]] <- add_data[[ii]] %>%
      filter_extreme(
        order = c(vars(temp_date), sources[[ii]]$order),
        by_vars = subject_keys,
        mode = sources[[ii]]$mode
      ) %>%
      mutate(
        temp_source_nr = ii,
        temp_date = !!sources[[ii]]$date,
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
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(
          !!!subject_keys,
          temp_source_nr,
          temp_date,
          DTHCAUS,
          !!!sources[[ii]]$traceability
        )
    } else {
      add_data[[ii]] <- add_data[[ii]] %>%
        select(!!!subject_keys, temp_source_nr, temp_date, DTHCAUS)
    }
  }
  # if a subject has multiple death info, keep the one from the first source
  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(
      order = vars(temp_date, temp_source_nr),
      by_vars = subject_keys,
      mode = "first"
    ) %>%
    select(-starts_with("temp_"))

  derive_vars_merged(dataset, dataset_add = dataset_add, by_vars = subject_keys)
}

#' Create a `dthcaus_source` Object
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the death cause.
#'
#' @param filter An expression used for filtering `dataset`.
#'
#' @param date A character vector to be used for sorting `dataset`.
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing.
#'
#'   If `NULL`: no date imputation is performed and partial dates are returned as missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'mm-dd': e.g. '06-15' for the 15th
#'   of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month.
#'
#'   Default is `NULL`
#'
#' @param order Sort order
#'
#'   Additional character vector to be used for sorting `dataset` to avoid duplicate record warning.
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
#' @param dataset Deprecated, please use `dataset_name` instead.
#'
#' @author Shimeng Huang
#'
#' @keywords source_specifications
#'
#' @seealso [derive_var_dthcaus()]
#'
#' @export
#'
#' @return An object of class "dthcaus_source".
dthcaus_source <- function(dataset_name,
                           filter,
                           date,
                           date_imputation = NULL,
                           order = NULL,
                           mode = "first",
                           dthcaus,
                           traceability_vars = NULL,
                           dataset = deprecated()) {
  if (!missing(dataset)) {
    deprecate_warn("0.6.0", "dthcaus_source(dataset = )", "dthcaus_source(dataset_name = )")
    dataset_name <- deparse(substitute(dataset))
  }

  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    date_imputation = assert_character_scalar(date_imputation, optional = TRUE),
    order = assert_order_vars(order, optional = TRUE),
    mode = assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE),
    dthcaus = assert_symbol(enquo(dthcaus)) %or% assert_character_scalar(dthcaus),
    traceability = assert_varval_list(traceability_vars, optional = TRUE)
  )
  class(out) <- c("dthcaus_source", "list")
  out
}
