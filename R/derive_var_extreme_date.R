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
#'   element are selected and observations where `date` is `NA` are removed.
#'   Then for each patient the first or last observation (with respect to `date`
#'   and `mode`) is selected.
#'
#'   1. The new variable is set to the variable specified by the `date` element.
#'   If this is a date variable (rather than datetime), then the time is imputed
#'   as `"00:00:00"`.
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
#' @family der_adsl
#' @keywords der_adsl
#'
#' @seealso [date_source()], [derive_var_extreme_dt()],
#'   [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data("admiral_dm")
#' data("admiral_ae")
#' data("admiral_lb")
#' data("admiral_adsl")
#'
#' # derive last known alive datetime (LSTALVDTM)
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = AESTDTM
#' )
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = AEENDTM
#' )
#'
#' ae_ext <- admiral_ae %>%
#'   derive_vars_dtm(
#'     dtc = AESTDTC,
#'     new_vars_prefix = "AEST",
#'     highest_imputation = "M"
#'   ) %>%
#'   derive_vars_dtm(
#'     dtc = AEENDTC,
#'     new_vars_prefix = "AEEN",
#'     highest_imputation = "M"
#'   )
#'
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = LBDTM,
#'   filter = !is.na(LBDTM)
#' )
#'
#' lb_ext <- derive_vars_dtm(
#'   admiral_lb,
#'   dtc = LBDTC,
#'   new_vars_prefix = "LB"
#' )
#'
#' adsl_date <- date_source(dataset_name = "adsl", date = TRTEDTM)
#'
#' admiral_dm %>%
#'   derive_var_extreme_dtm(
#'     new_var = LSTALVDTM,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = admiral_adsl,
#'       ae = ae_ext, lb = lb_ext
#'     ),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDTM)
#'
#' # derive last alive datetime and traceability variables
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = AESTDTM,
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = AEENDTM,
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = LBDTM,
#'   filter = !is.na(LBDTM),
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
#'     source_datasets = list(
#'       adsl = admiral_adsl,
#'       ae = ae_ext,
#'       lb = lb_ext
#'     ),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDTM, LALVDOM, LALVSEQ, LALVVAR)
derive_var_extreme_dtm <- function(dataset,
                                   new_var,
                                   ...,
                                   source_datasets,
                                   mode,
                                   subject_keys = get_admiral_option("subject_keys")) {
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
    assert_date_var(
      dataset = source_dataset,
      var = !!date,
      dataset_name = source_dataset_name
    )

    if (!is.null(sources[[i]]$traceability_vars)) {
      warn_if_vars_exist(source_dataset, names(sources[[i]]$traceability_vars))
      assert_data_frame(
        source_dataset,
        required_vars = get_source_vars(sources[[i]]$traceability_vars)
      )
    }

    add_data[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      filter(!is.na(!!date)) %>%
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
      !!new_var := convert_date_to_dtm(!!date)
    )
  }

  all_data <- add_data %>%
    bind_rows() %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(!!new_var),
      mode = mode,
      check_type = "none"
    )

  derive_vars_merged(
    dataset,
    dataset_add = all_data,
    by_vars = subject_keys
  )
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
#'   element are selected and observations where `date` is `NA` are removed.
#'   Then for each patient the first or last observation (with respect to `date`
#'   and `mode`) is selected.
#'
#'   1. The new variable is set to the variable specified by the `date` element.
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
#' @family der_adsl
#' @keywords der_adsl
#'
#' @seealso [date_source()], [derive_var_extreme_dtm()], [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data("admiral_dm")
#' data("admiral_ae")
#' data("admiral_lb")
#' data("admiral_adsl")
#'
#' # derive last known alive date (LSTALVDT)
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = AESTDT
#' )
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = AEENDT
#' )
#'
#' ae_ext <- admiral_ae %>%
#'   derive_vars_dt(
#'     dtc = AESTDTC,
#'     new_vars_prefix = "AEST",
#'     highest_imputation = "M"
#'   ) %>%
#'   derive_vars_dt(
#'     dtc = AEENDTC,
#'     new_vars_prefix = "AEEN",
#'     highest_imputation = "M"
#'   )
#'
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = LBDT,
#'   filter = !is.na(LBDT),
#' )
#'
#' lb_ext <- derive_vars_dt(
#'   admiral_lb,
#'   dtc = LBDTC,
#'   new_vars_prefix = "LB"
#' )
#'
#' adsl_date <- date_source(dataset_name = "adsl", date = TRTEDT)
#'
#' admiral_dm %>%
#'   derive_var_extreme_dt(
#'     new_var = LSTALVDT,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = admiral_adsl,
#'       ae = ae_ext,
#'       lb = lb_ext
#'     ),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDT)
#'
#' # derive last alive date and traceability variables
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = AESTDT,
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = AEENDT,
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = LBDT,
#'   filter = !is.na(LBDT),
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
#'     source_datasets = list(
#'       adsl = admiral_adsl,
#'       ae = ae_ext,
#'       lb = lb_ext
#'     ),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR)
derive_var_extreme_dt <- function(dataset,
                                  new_var,
                                  ...,
                                  source_datasets,
                                  mode,
                                  subject_keys = get_admiral_option("subject_keys")) {
  new_var <- assert_symbol(enquo(new_var))

  sources <- list(...)
  assert_list_of(sources, "date_source")

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
#' @param date A variable providing a date. A date or a datetime can be
#'   specified. An unquoted symbol is expected.
#'
#' @param date_imputation *Deprecated*, please use `derive_vars_dtm()` to
#'   convert DTC variables to datetime variables in the dataset.
#'
#' @param time_imputation *Deprecated*, please use `derive_vars_dtm()` to
#'   convert DTC variables to datetime variables in the dataset.
#'
#' @param preserve *Deprecated*, please use `derive_vars_dtm()` to convert DTC
#'   variables to datetime variables in the dataset.
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
#' @family source_specifications
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class `date_source`.
#'
#' @examples
#'
#' # treatment end date from ADSL
#' trt_end_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDT
#' )
#'
#' # lab date from LB where assessment was taken, i.e. not "NOT DONE"
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   filter = LBSTAT != "NOT DONE" | is.na(LBSTAT),
#'   date = LBDT
#' )
#'
#' # death date from ADSL including traceability variables
#' death_date <- date_source(
#'   dataset_name = "adsl",
#'   date = DTHDT,
#'   traceability_vars = vars(
#'     LALVDOM = "ADSL",
#'     LALVVAR = "DTHDT"
#'   )
#' )
date_source <- function(dataset_name,
                        filter = NULL,
                        date,
                        date_imputation = deprecated(),
                        time_imputation = deprecated(),
                        preserve = deprecated(),
                        traceability_vars = NULL) {
  if (!missing(date_imputation)) {
    deprecate_stop(
      "0.8.0",
      "date_source(date_imputation = )",
      details = paste0(
        "Please use `derive_vars_dtm()` to convert DTC variables",
        " to datetime variables in the dataset."
      )
    )
  }
  if (!missing(time_imputation)) {
    deprecate_stop(
      "0.8.0",
      "date_source(time_imputation = )",
      details = paste0(
        "Please use `derive_vars_dtm()` to convert DTC variables",
        " to datetime variables in the dataset."
      )
    )
  }
  if (!missing(preserve)) {
    deprecate_stop(
      "0.8.0",
      "date_source(preserve = )",
      details = paste0(
        "Please use `derive_vars_dtm()` to convert DTC variables",
        " to datetime variables in the dataset."
      )
    )
  }
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    traceability_vars = assert_varval_list(traceability_vars, optional = TRUE)
  )
  class(out) <- c("date_source", "source", "list")
  out
}
