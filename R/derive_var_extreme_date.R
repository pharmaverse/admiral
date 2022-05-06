#' Derive Last Known Alive Date
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' *Deprecated*, please use `derive_var_extreme_dt()` instead. Add the last
#' known alive date (`LSTALVDT`) to the dataset.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by `subject_keys` are required.
#'
#' @param source_datasets A named `list` containing datasets in which to search for the
#'   last known alive date
#'
#' @param ... Source(s) of known alive dates. One or more `lstalvdt_source()` objects are
#'   expected.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @author Stefan Bundfuss, Thomas Neitmann
#'
#' @return The input dataset with the `LSTALVDT` variable added.
#'
#' @keywords derivation adsl
#'
#' @export
derive_var_lstalvdt <- function(dataset,
                                ...,
                                source_datasets,
                                subject_keys = vars(STUDYID, USUBJID)) {
  deprecate_warn("0.7.0", "derive_var_lstalvdt()", "derive_var_extreme_dt()")
  derive_var_extreme_dt(
    dataset,
    new_var = LSTALVDT,
    ...,
    source_datasets = source_datasets,
    mode = "last",
    subject_keys = subject_keys
  )
}

#' Derive First or Last Datetime from Multiple Sources
#'
#' Add the first or last datetime from multiple sources to the dataset, e.g.,
#' the last known alive datetime (`LSTALVDTM`).
#'
#' @param dataset Input dataset
#'
#'   The variables specified by `subject_keys` are required.
#'
#' @param new_var Name of variable to create
#'
#' @param source_datasets A named `list` containing datasets in which to search
#'   for the first or last date
#'
#' @param ... Source(s) of dates. One or more `date_source()` objects are
#'   expected.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first date for each subject is selected. If
#'   `"last"` is specified, the last date for each subject is selected.
#'
#'   Permitted Values:  `"first"`, `"last"`
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   1. For each source dataset the observations as specified by the `filter`
#'   element are selected. Then for each patient the first or last observation
#'   (with respect to `date` and `mode`) is selected.
#'
#'   1. The new variable is set to the variable specified by the `date` element.
#'   If the date variable is a date variable, the time is imputed as specified
#'   by the `time_imputation` element. If the source variable is a character
#'   variable, it is converted to a datetime. If the date is incomplete, it is
#'   imputed as specified by the `date_imputation` and `time_imputation`
#'   element.
#'
#'   1. The variables specified by the `traceability_vars` element are added.
#'
#'   1. The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   1. For each patient the first or last observation (with respect to the new
#'   variable and `mode`) from the single dataset is selected and the new
#'   variable is merged to the input dataset.
#'
#' @return The input dataset with the new variable added.
#'
#' @author Stefan Bundfuss, Thomas Neitmann
#'
#' @keywords derivation adsl
#'
#' @seealso [date_source()], [derive_var_extreme_dt()],
#'   [derive_vars_merged_dt()], [derive_vars_merged_dtm()],
#'   [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiraltest)
#' data("admiral_dm")
#' data("admiral_ae")
#' data("admiral_lb")
#' data("adsl")
#'
#' # derive last known alive datetime (LSTALVDTM)
#' ae_start <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AESTDTC,
#'   date_imputation = "first",
#'   time_imputation = "first"
#' )
#' ae_end <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AEENDTC,
#'   date_imputation = "first",
#'   time_imputation = "first"
#' )
#' lb_date <- date_source(
#'   dataset_name = "admiral_lb",
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10,
#'   time_imputation = "first"
#' )
#' adsl_date <- date_source(dataset_name = "adsl", date = TRTEDTM)
#'
#' admiral_dm %>%
#'   derive_var_extreme_dtm(
#'     new_var = LSTALVDTM,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(adsl = adsl, admiral_ae = admiral_ae, admiral_lb = admiral_lb),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDTM)
#'
#' # derive last alive datetime and traceability variables
#' ae_start <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AESTDTC,
#'   date_imputation = "first",
#'   time_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AEENDTC,
#'   date_imputation = "first",
#'   time_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- date_source(
#'   dataset_name = "admiral_lb",
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10,
#'   time_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDTM,
#'   traceability_vars = vars(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDTM"
#'   )
#' )
#'
#' admiral_dm %>%
#'   derive_var_extreme_dtm(
#'     new_var = LSTALVDTM,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(adsl = adsl, admiral_ae = admiral_ae, admiral_lb = admiral_lb),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDTM, LALVDOM, LALVSEQ, LALVVAR)
derive_var_extreme_dtm <- function(dataset,
                                   new_var,
                                   ...,
                                   source_datasets,
                                   mode,
                                   subject_keys = vars(STUDYID, USUBJID)) {
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)
  new_var <- assert_symbol(enquo(new_var))
  assert_list_of(source_datasets, "data.frame")
  sources <- rlang::list2(...)
  assert_list_of(sources, "date_source")
  mode <- assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE
  )

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

  warn_if_vars_exist(dataset, vars2chr(quo_c(new_var)))

  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (i > 1) {
      warn_if_inconsistent_list(
        base = sources[[i - 1]]$traceability_vars,
        compare = sources[[i]]$traceability_vars,
        list_name = "date_source()",
        i = i
      )
    }

    source_dataset_name <- sources[[i]]$dataset_name
    source_dataset <- source_datasets[[source_dataset_name]]

    date <- quo_get_expr(sources[[i]]$date)
    add_data[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      filter_extreme(
        order = vars(!!date),
        by_vars = subject_keys,
        mode = mode,
        check_type = "none"
      )

    add_data[[i]] <- transmute(
      add_data[[i]],
      !!!subject_keys,
      !!!sources[[i]]$traceability_vars,
      !!new_var := convert_date_to_dtm(
        !!date,
        date_imputation = sources[[i]]$date_imputation,
        time_imputation = sources[[i]]$time_imputation,
        preserve = sources[[i]]$preserve
      )
    )
  }

  all_data <- add_data %>%
    bind_rows() %>%
    filter(!is.na(!!new_var)) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(!!new_var),
      mode = mode,
      check_type = "none"
    )

  derive_vars_merged(dataset,
                     dataset_add = all_data,
                     by_vars = subject_keys)
}

#' Derive First or Last Date from Multiple Sources
#'
#' Add the first or last date from multiple sources to the dataset, e.g.,
#' the last known alive date (`LSTALVDT`).
#'
#' @inheritParams derive_var_extreme_dtm
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   1. For each source dataset the observations as specified by the `filter`
#'   element are selected. Then for each patient the first or last observation
#'   (with respect to `date` and `mode`) is selected.
#'
#'   1. The new variable is set to the variable specified by the `date` element.
#'   If the date variable is a date variable, the time is imputed as
#'   `time_imputation = "first"`. If the source variable is a character
#'   variable, it is converted to a datetime. If the date is incomplete, it is
#'   imputed as specified by the `date_imputation` element and with
#'   `time_imputation = "first"`.
#'
#'   1. The variables specified by the `traceability_vars` element are added.
#'
#'   1. The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   1. For each patient the first or last observation (with respect to the new
#'   variable and `mode`) from the single dataset is selected and the new
#'   variable is merged to the input dataset.
#'
#'   1. The time part is removed from the new variable.
#'
#' @return The input dataset with the new variable added.
#'
#' @author Stefan Bundfuss, Thomas Neitmann
#'
#' @keywords derivation adsl
#'
#' @seealso [date_source()], [derive_var_extreme_dtm()],
#'   [derive_vars_merged_dt()], [derive_vars_merged_dtm()],
#'   [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiraltest)
#' data("admiral_dm")
#' data("admiral_ae")
#' data("admiral_lb")
#' data("adsl")
#'
#' # derive last known alive date (LSTALVDT)
#' ae_start <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AESTDTC,
#'   date_imputation = "first",
#' )
#' ae_end <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AEENDTC,
#'   date_imputation = "first",
#' )
#' lb_date <- date_source(
#'   dataset_name = "admiral_lb",
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10,
#' )
#' adsl_date <- date_source(dataset_name = "adsl", date = TRTEDT)
#'
#' admiral_dm %>%
#'   derive_var_extreme_dt(
#'     new_var = LSTALVDT,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(adsl = adsl, admiral_ae = admiral_ae, admiral_lb = admiral_lb),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDT)
#'
#' # derive last alive date and traceability variables
#' ae_start <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AESTDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "admiral_ae",
#'   date = AEENDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- date_source(
#'   dataset_name = "admiral_lb",
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10,
#'   traceability_vars = vars(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDT,
#'   traceability_vars = vars(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDT"
#'   )
#' )
#'
#' admiral_dm %>%
#'   derive_var_extreme_dt(
#'     new_var = LSTALVDT,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(adsl = adsl, admiral_ae = admiral_ae, admiral_lb = admiral_lb),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR)
derive_var_extreme_dt <- function(dataset,
                                  new_var,
                                  ...,
                                  source_datasets,
                                  mode,
                                  subject_keys = vars(STUDYID, USUBJID)) {
  new_var <- assert_symbol(enquo(new_var))

  sources <- list(...)
  assert_list_of(sources, "date_source")
  for (i in seq_along(sources)) {
    sources[[i]]$time_imputation = "first"
  }

  derive_var_extreme_dtm(
    dataset,
    new_var = !!new_var,
    !!!sources,
    source_datasets = source_datasets,
    mode = mode,
    subject_keys = subject_keys
  ) %>%
    mutate(!!new_var := date(!!new_var))
}

#' Create an `lstalvdt_source` object
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' *Deprecated*, please use `date_source()` instead.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the last known alive date.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date A variable providing a date where the patient was known to be
#'   alive. A date, a datetime, or a character variable containing ISO 8601
#'   dates can be specified. An unquoted symbol is expected.
#'
#' @param date_imputation A string defining the date imputation for `date`.
#'   See `date_imputation` parameter of `impute_dtc()` for valid values.
#'
#' @param traceability_vars A named list returned by `vars()` defining the
#'   traceability variables, e.g. `vars(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR
#'   = "AESTDTC")`. The values must be a symbol, a character string, or `NA`.
#'
#' @param dataset Deprecated, please use `dataset_name` instead.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class `lstalvdt_source`.
lstalvdt_source <- function(dataset_name,
                            filter = NULL,
                            date,
                            date_imputation = NULL,
                            traceability_vars = NULL,
                            dataset = deprecated()) {
  if (!missing(dataset)) {
    deprecate_warn("0.6.0", "lstalvdt_source(dataset = )", "lstalvdt_source(dataset_name = )")
    dataset_name <- deparse(substitute(dataset))
  }
  deprecate_warn("0.7.0", "lstalvdt_source()", "date_source()")

  date_source(
    dataset_name = dataset_name,
    filter = !!enquo(filter),
    date = !!enquo(date),
    date_imputation = date_imputation,
    traceability_vars = traceability_vars)
}

#' Create a `date_source` object
#'
#' Create a `date_source` object as input for `derive_var_extreme_dt()` and
#' `derive_var_extreme_dtm()`.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the date.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date A variable providing a date. A date, a datetime, or a character
#'   variable containing ISO 8601 dates can be specified. An unquoted symbol is
#'   expected.
#'
#' @param date_imputation A string defining the date imputation for `date`.
#'   See `date_imputation` parameter of `impute_dtc()` for valid values.
#'
#' @param time_imputation A string defining the time imputation for `date`.
#'   See `time_imputation` parameter of `impute_dtc()` for valid values.
#'
#' @param preserve Should day be preserved if month is imputed for `date`.
#'   See `preserve` parameter of `impute_dtc()` for details.
#'
#' @param traceability_vars A named list returned by `vars()` defining the
#'   traceability variables, e.g. `vars(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR
#'   = "AESTDTC")`. The values must be a symbol, a character string, a numeric,
#'   or `NA`.
#'
#' @author Stefan Bundfuss
#'
#' @seealso [derive_var_extreme_dtm()], [derive_var_extreme_dt()]
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class `date_source`.
date_source <- function(dataset_name,
                        filter = NULL,
                        date,
                        date_imputation = NULL,
                        time_imputation = NULL,
                        preserve = FALSE,
                        traceability_vars = NULL) {

  if (!is.null(date_imputation)) {
    assert_that(is_valid_date_entry(date_imputation))
  }
  if (!is.null(time_imputation)) {
    assert_that(is_valid_time_entry(time_imputation))
  }
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    date_imputation = date_imputation,
    time_imputation = time_imputation,
    preserve = preserve,
    traceability_vars = assert_varval_list(traceability_vars, optional = TRUE)
  )
  class(out) <- c("date_source", "list")
  out
}
