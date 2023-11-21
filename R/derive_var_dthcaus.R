#' Derive Death Cause
#'
#' @description
#' `r lifecycle::badge("superseded")` The `derive_var_dthcaus()`
#' function has been superseded in favor of `derive_vars_extreme_event()`.
#'
#' Derive death cause (`DTHCAUS`) and add traceability variables if required.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("subject_keys"))`
#'
#' @param source_datasets A named `list` containing datasets in which to search for the
#'   death cause
#'
#' @param ... Objects of class "dthcaus_source" created by [`dthcaus_source()`].
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of expressions where the expressions are symbols as returned by
#' `exprs()` is expected.
#'
#' @details
#' This function derives `DTHCAUS` along with the user-defined traceability
#' variables, if required. If a subject has death info from multiple sources,
#' the one from the source with the earliest death date will be used. If dates are
#' equivalent, the first source will be kept, so the user should provide the inputs in
#' the preferred order.
#'
#' @family superseded
#' @keywords superseded
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
#' )
#'
#' ds <- tribble(
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
#'   date = convert_dtc_to_dt(AEDTHDTC),
#'   mode = "first",
#'   dthcaus = AEDECOD
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = convert_dtc_to_dt(DSSTDTC),
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
#'   date = convert_dtc_to_dt(AEDTHDTC),
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = convert_dtc_to_dt(DSSTDTC),
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(adsl, src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))
#'
#' # Derive `DTHCAUS` as above - now including post-study deaths with different `DTHCAUS` value
#' src_ae <- dthcaus_source(
#'   dataset_name = "ae",
#'   filter = AEOUT == "FATAL",
#'   date = convert_dtc_to_dt(AEDTHDTC),
#'   mode = "first",
#'   dthcaus = AEDECOD,
#'   set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
#' )
#'
#' ds <- mutate(
#'   ds,
#'   DSSTDT = convert_dtc_to_dt(DSSTDTC)
#' )
#'
#' src_ds <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = DSTERM,
#'   set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' src_ds_post <- dthcaus_source(
#'   dataset_name = "ds",
#'   filter = DSDECOD == "DEATH" & DSTERM == "POST STUDY REPORTING OF DEATH",
#'   date = DSSTDT,
#'   mode = "first",
#'   dthcaus = "POST STUDY: UNKNOWN CAUSE",
#'   set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
#' )
#'
#' derive_var_dthcaus(
#'   adsl,
#'   src_ae, src_ds, src_ds_post,
#'   source_datasets = list(ae = ae, ds = ds)
#' )
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
  tmp_source_nr <- get_new_tmp_var(dataset)
  tmp_date <- get_new_tmp_var(dataset)
  for (ii in seq_along(sources)) {
    source_dataset_name <- sources[[ii]]$dataset_name
    source_dataset <- source_datasets[[source_dataset_name]]
    add_data[[ii]] <- source_dataset %>%
      filter_if(sources[[ii]]$filter)

    date <- sources[[ii]]$date
    if (is.symbol(date)) {
      date_var <- date
    } else {
      date_var <- get_new_tmp_var(dataset = add_data[[ii]], prefix = "tmp_date")
      add_data[[ii]] <- mutate(
        add_data[[ii]],
        !!date_var := !!date
      )
    }
    assert_date_var(
      dataset = add_data[[ii]],
      var = !!date_var,
      dataset_name = source_dataset_name
    )

    # if several death records, use the first/last according to 'mode'
    add_data[[ii]] <- add_data[[ii]] %>%
      filter_extreme(
        order = exprs(!!date_var, !!!sources[[ii]]$order),
        by_vars = subject_keys,
        mode = sources[[ii]]$mode
      ) %>%
      mutate(
        !!tmp_source_nr := ii,
        !!tmp_date := !!date_var,
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
      order = exprs(!!tmp_date, !!tmp_source_nr),
      by_vars = subject_keys,
      mode = "first"
    ) %>%
    remove_tmp_vars()

  derive_vars_merged(dataset, dataset_add = dataset_add, by_vars = subject_keys)
}

#' Create a `dthcaus_source` Object
#'
#' @description
#' `r lifecycle::badge("superseded")` The `derive_var_dthcaus()`
#' function and `dthcaus_source()` have been superseded in favor of
#' `derive_vars_extreme_event()`.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the death cause.
#'
#' @param filter An expression used for filtering `dataset`.
#'
#' @param date A date or datetime variable or an expression to be used for
#'   sorting `dataset`.
#'
#' @param order Sort order
#'
#'   Additional variables/expressions to be used for sorting the `dataset`. The
#'   dataset is ordered by `date` and `order`. Can be used to avoid duplicate
#'   record warning.
#'
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @param mode One of `"first"` or `"last"`.
#' Either the `"first"` or `"last"` observation is preserved from the `dataset`
#' which is ordered by `date`.
#'
#' @param dthcaus A variable name, an expression, or a string literal
#'
#'   If a variable name is specified, e.g., `AEDECOD`, it is the variable in the
#'   source dataset to be used to assign values to `DTHCAUS`; if an expression,
#'   e.g., `str_to_upper(AEDECOD)`, it is evaluated in the source dataset and
#'   the results is assigned to `DTHCAUS`; if a string literal, e.g. `"Adverse
#'   Event"`, it is the fixed value to be assigned to `DTHCAUS`.
#'
#'
#' @param traceability_vars A named list returned by [`exprs()`] listing the
#'   traceability variables, e.g. `exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)`. The
#'   left-hand side (names of the list elements) gives the names of the
#'   traceability variables in the returned dataset. The right-hand side (values
#'   of the list elements) gives the values of the traceability variables in the
#'   returned dataset. These can be either strings, numbers, symbols, or
#'   expressions referring to existing variables.
#'
#'    `r lifecycle::badge("deprecated")` Please use `set_values_to` instead.
#'
#' @param set_values_to Variables to be set to trace the source dataset
#'
#' @family superseded
#' @keywords superseded
#'
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
#'   date = convert_dtc_to_dt(DSSTDTC),
#'   mode = "first",
#'   dthcaus = DSTERM
#' )
dthcaus_source <- function(dataset_name,
                           filter,
                           date,
                           order = NULL,
                           mode = "first",
                           dthcaus,
                           set_values_to = NULL,
                           traceability_vars = NULL) {
  if (!is.null(traceability_vars)) {
    deprecate_stop(
      "0.12.0",
      "dthcaus_source(traceability_vars = )",
      "dthcaus_source(set_values_to = )"
    )
    set_values_to <- traceability_vars
  }

  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enexpr(filter), optional = TRUE),
    date = assert_expr(enexpr(date)),
    order = assert_expr_list(order, optional = TRUE),
    mode = assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE),
    dthcaus = assert_expr(enexpr(dthcaus)),
    traceability = assert_expr_list(set_values_to, named = TRUE, optional = TRUE)
  )
  class(out) <- c("dthcaus_source", "source", "list")
  out
}
