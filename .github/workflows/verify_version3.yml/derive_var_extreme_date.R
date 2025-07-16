#' Derive First or Last Datetime from Multiple Sources
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_var_extreme_dtm()`
#' function has been deprecated in favor of `derive_vars_extreme_event()`.
#'
#' Add the first or last datetime from multiple sources to the dataset, e.g.,
#' the last known alive datetime (`LSTALVDTM`).
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("subject_keys"))`
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
#' @permitted  `"first"`, `"last"`
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of expressions where the expressions are symbols as returned by
#'   `exprs()` is expected.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   1. For each source dataset the observations as specified by the `filter`
#'   element are selected and observations where `date` is `NA` are removed.
#'   Then for each patient the first or last observation (with respect to `date`
#'   and `mode`) is selected.
#'
#'   1. The new variable is set to the variable or expression specified by the
#'   `date` element. If this is a date variable (rather than datetime), then the
#'   time is imputed as `"00:00:00"`.
#'
#'   1. The variables specified by the `set_values_to` element are added.
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
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [date_source()], [derive_var_extreme_dt()],
#'   [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#' dm <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",    "DM", "01-1130",   84, "YEARS",
#'   "PILOT01",    "DM", "01-1133",   81, "YEARS",
#'   "PILOT01",    "DM", "01-1211",   76, "YEARS",
#'   "PILOT01",    "DM", "09-1081",   86, "YEARS",
#'   "PILOT01",    "DM", "09-1088",   69, "YEARS"
#' )
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AESEQ,     ~AESTDTC,     ~AEENDTC,
#'   "PILOT01",    "AE", "01-1130",      5, "2014-05-09", "2014-05-09",
#'   "PILOT01",    "AE", "01-1130",      6, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      4, "2014-05-09", "2014-05-09",
#'   "PILOT01",    "AE", "01-1130",      8, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      7, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      2, "2014-03-09", "2014-03-09",
#'   "PILOT01",    "AE", "01-1130",      1, "2014-03-09", "2014-03-16",
#'   "PILOT01",    "AE", "01-1130",      3, "2014-03-09", "2014-03-16",
#'   "PILOT01",    "AE", "01-1133",      1, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      3, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      2, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      4, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1211",      5, "2012-11-29",           NA,
#'   "PILOT01",    "AE", "01-1211",      1, "2012-11-16",           NA,
#'   "PILOT01",    "AE", "01-1211",      7, "2013-01-11",           NA,
#'   "PILOT01",    "AE", "01-1211",      8, "2013-01-11",           NA,
#'   "PILOT01",    "AE", "01-1211",      4, "2012-11-22",           NA,
#'   "PILOT01",    "AE", "01-1211",      2, "2012-11-21", "2012-11-21",
#'   "PILOT01",    "AE", "01-1211",      3, "2012-11-21",           NA,
#'   "PILOT01",    "AE", "01-1211",      6, "2012-12-09",           NA,
#'   "PILOT01",    "AE", "01-1211",      9, "2013-01-14", "2013-01-14",
#'   "PILOT01",    "AE", "09-1081",      2, "2014-05-01",           NA,
#'   "PILOT01",    "AE", "09-1081",      1, "2014-04-07",           NA,
#'   "PILOT01",    "AE", "09-1088",      1, "2014-05-08",           NA,
#'   "PILOT01",    "AE", "09-1088",      2, "2014-08-02",           NA
#' )
#' lb <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
#'   "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
#'   "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
#'   "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
#'   "PILOT01",    "LB", "01-1133",    304, "2013-04-29T10:13",
#'   "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
#'   "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
#'   "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
#'   "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
#'   "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
#'   "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
#' )
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID,              ~TRTEDTM,
#'   "PILOT01", "01-1130", "2014-08-16 23:59:59",
#'   "PILOT01", "01-1133", "2013-04-28 23:59:59",
#'   "PILOT01", "01-1211", "2013-01-12 23:59:59",
#'   "PILOT01", "09-1081", "2014-04-27 23:59:59",
#'   "PILOT01", "09-1088", "2014-10-09 23:59:59"
#' ) %>%
#'   mutate(
#'     TRTEDTM = as_datetime(TRTEDTM)
#'   )
#'
#' # derive last known alive datetime (LSTALVDTM)
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dtm(AESTDTC, highest_imputation = "M"),
#' )
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dtm(AEENDTC, highest_imputation = "M"),
#' )
#'
#' ae_ext <- ae %>%
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
#'   date = convert_dtc_to_dtm(LBDTC),
#' )
#'
#' lb_ext <- derive_vars_dtm(
#'   lb,
#'   dtc = LBDTC,
#'   new_vars_prefix = "LB"
#' )
#'
#' adsl_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDTM
#' )
#'
#' dm %>%
#'   derive_var_extreme_dtm(
#'     new_var = LSTALVDTM,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = adsl,
#'       ae = ae_ext,
#'       lb = lb_ext
#'     ),
#'     mode = "last"
#'   ) %>%
#'   select(USUBJID, LSTALVDTM)
#'
#' # derive last alive datetime and traceability variables
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dtm(AESTDTC, highest_imputation = "M"),
#'   set_values_to = exprs(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dtm(AEENDTC, highest_imputation = "M"),
#'   set_values_to = exprs(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = convert_dtc_to_dtm(LBDTC),
#'   set_values_to = exprs(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDTM,
#'   set_values_to = exprs(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDTM"
#'   )
#' )
#'
#' dm %>%
#'   derive_var_extreme_dtm(
#'     new_var = LSTALVDTM,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = adsl,
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
  deprecate_inform(
    when = "1.2.0",
    what = "derive_var_extreme_dtm()",
    with = "derive_vars_extreme_event()",
    details = c(
      x = "This message will turn into a warning at the beginning of 2026.",
      i = "See admiral's deprecation guidance:
      https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )

  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)
  new_var <- assert_symbol(enexpr(new_var))
  assert_list_of(source_datasets, "data.frame")
  sources <- list2(...)
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
    message_text = c(
      paste0(
        "The dataset names must be included in the list specified for the ",
        "{.arg source_datasets} argument."
      ),
      i = paste(
        "Following names were provided by {.arg source_datasets}:",
        ansi_collapse(source_names)
      )
    )
  )

  warn_if_vars_exist(dataset, vars2chr(expr_c(new_var)))

  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (i > 1) {
      warn_if_inconsistent_list(
        base = sources[[i - 1]]$set_values_to,
        compare = sources[[i]]$set_values_to,
        list_name = "date_source()",
        i = i
      )
    }

    source_dataset_name <- sources[[i]]$dataset_name
    source_dataset <- source_datasets[[source_dataset_name]]

    date <- sources[[i]]$date
    if (is.symbol(date)) {
      date_var <- date
    } else {
      date_var <- get_new_tmp_var(dataset = source_dataset, prefix = "tmp_date")
      source_dataset <- mutate(
        source_dataset,
        !!date_var := !!date
      )
    }
    assert_date_var(
      dataset = source_dataset,
      var = !!date_var,
      dataset_name = source_dataset_name
    )

    if (!is.null(sources[[i]]$set_values_to)) {
      warn_if_vars_exist(source_dataset, names(sources[[i]]$set_values_to))
      assert_data_frame(
        source_dataset,
        required_vars = get_source_vars(sources[[i]]$set_values_to)
      )
    }

    add_data[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      filter(!is.na(!!date_var)) %>%
      filter_extreme(
        order = exprs(!!date_var),
        by_vars = subject_keys,
        mode = mode,
        check_type = "none"
      )

    add_data[[i]] <- mutate(
      add_data[[i]],
      !!!subject_keys,
      !!!sources[[i]]$set_values_to,
      !!new_var := convert_date_to_dtm(!!date_var),
      .keep = "none"
    )
  }

  all_data <- add_data %>%
    bind_rows() %>%
    filter_extreme(
      by_vars = subject_keys,
      order = exprs(!!new_var),
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
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_var_extreme_dt()`
#' function has been deprecated in favor of `derive_vars_extreme_event()`.
#'
#' Add the first or last date from multiple sources to the
#' dataset, e.g., the last known alive date (`LSTALVDT`).
#'
#' **Note:** This is a wrapper function for the function `derive_var_extreme_dtm()`.
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
#'   1. The new variable is set to the variable or expression specified by the
#'   `date` element.
#'
#'   1. The variables specified by the `set_values_to` element are added.
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
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [date_source()], [derive_var_extreme_dtm()], [derive_vars_merged()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AESEQ,     ~AESTDTC,     ~AEENDTC,
#'   "PILOT01",    "AE", "01-1130",      5, "2014-05-09", "2014-05-09",
#'   "PILOT01",    "AE", "01-1130",      6, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      4, "2014-05-09", "2014-05-09",
#'   "PILOT01",    "AE", "01-1130",      8, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      7, "2014-05-22",           NA,
#'   "PILOT01",    "AE", "01-1130",      2, "2014-03-09", "2014-03-09",
#'   "PILOT01",    "AE", "01-1130",      1, "2014-03-09", "2014-03-16",
#'   "PILOT01",    "AE", "01-1130",      3, "2014-03-09", "2014-03-16",
#'   "PILOT01",    "AE", "01-1133",      1, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      3, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      2, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1133",      4, "2012-12-27",           NA,
#'   "PILOT01",    "AE", "01-1211",      5, "2012-11-29",           NA,
#'   "PILOT01",    "AE", "01-1211",      1, "2012-11-16",           NA,
#'   "PILOT01",    "AE", "01-1211",      7, "2013-01-11",           NA,
#'   "PILOT01",    "AE", "01-1211",      8, "2013-01-11",           NA,
#'   "PILOT01",    "AE", "01-1211",      4, "2012-11-22",           NA,
#'   "PILOT01",    "AE", "01-1211",      2, "2012-11-21", "2012-11-21",
#'   "PILOT01",    "AE", "01-1211",      3, "2012-11-21",           NA,
#'   "PILOT01",    "AE", "01-1211",      6, "2012-12-09",           NA,
#'   "PILOT01",    "AE", "01-1211",      9, "2013-01-14", "2013-01-14",
#'   "PILOT01",    "AE", "09-1081",      2, "2014-05-01",           NA,
#'   "PILOT01",    "AE", "09-1081",      1, "2014-04-07",           NA,
#'   "PILOT01",    "AE", "09-1088",      1, "2014-05-08",           NA,
#'   "PILOT01",    "AE", "09-1088",      2, "2014-08-02",           NA
#' )
#'
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID,              ~TRTEDTM,      ~TRTEDT,
#'   "PILOT01", "01-1130", "2014-08-16 23:59:59", "2014-08-16",
#'   "PILOT01", "01-1133", "2013-04-28 23:59:59", "2013-04-28",
#'   "PILOT01", "01-1211", "2013-01-12 23:59:59", "2013-01-12",
#'   "PILOT01", "09-1081", "2014-04-27 23:59:59", "2014-04-27",
#'   "PILOT01", "09-1088", "2014-10-09 23:59:59", "2014-10-09"
#' ) %>%
#'   mutate(
#'     across(TRTEDTM:TRTEDT, as.Date)
#'   )
#'
#'
#' lb <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
#'   "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
#'   "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
#'   "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
#'   "PILOT01",    "LB", "01-1133",    304, "2013-04-29T10:13",
#'   "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
#'   "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
#'   "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
#'   "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
#'   "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
#'   "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
#' )
#'
#' dm <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",    "DM", "01-1130",   84, "YEARS",
#'   "PILOT01",    "DM", "01-1133",   81, "YEARS",
#'   "PILOT01",    "DM", "01-1211",   76, "YEARS",
#'   "PILOT01",    "DM", "09-1081",   86, "YEARS",
#'   "PILOT01",    "DM", "09-1088",   69, "YEARS"
#' )
#'
#' ae_start <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dt(AESTDTC, highest_imputation = "M")
#' )
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dt(AEENDTC, highest_imputation = "M")
#' )
#'
#' ae_ext <- ae %>%
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
#'   date = convert_dtc_to_dt(LBDTC)
#' )
#'
#' lb_ext <- derive_vars_dt(
#'   lb,
#'   dtc = LBDTC,
#'   new_vars_prefix = "LB"
#' )
#'
#' adsl_date <- date_source(dataset_name = "adsl", date = TRTEDT)
#'
#' dm %>%
#'   derive_var_extreme_dt(
#'     new_var = LSTALVDT,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = adsl,
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
#'   date = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
#'   set_values_to = exprs(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- date_source(
#'   dataset_name = "ae",
#'   date = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
#'   set_values_to = exprs(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#'
#' lb_date <- date_source(
#'   dataset_name = "lb",
#'   date = convert_dtc_to_dt(LBDTC),
#'   set_values_to = exprs(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- date_source(
#'   dataset_name = "adsl",
#'   date = TRTEDT,
#'   set_values_to = exprs(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDT"
#'   )
#' )
#'
#' dm %>%
#'   derive_var_extreme_dt(
#'     new_var = LSTALVDT,
#'     ae_start, ae_end, lb_date, adsl_date,
#'     source_datasets = list(
#'       adsl = adsl,
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
  deprecate_inform(
    when = "1.2.0",
    what = "derive_var_extreme_dt()",
    with = "derive_vars_extreme_event()",
    details = c(
      x = "This message will turn into a warning at the beginning of 2026.",
      i = "See admiral's deprecation guidance:
      https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )


  new_var <- assert_symbol(enexpr(new_var))

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
#' @description
#' `r lifecycle::badge("deprecated")` The `date_source()`
#' function has been deprecated in favor of `event()`.
#'
#' Create a `date_source` object as input for `derive_var_extreme_dt()` and
#' `derive_var_extreme_dtm()`.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the date.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date A variable or an expression providing a date. A date or a
#'   datetime can be specified. An unquoted symbol or expression is expected.
#'
#' @param set_values_to Variables to be set
#'
#' @seealso [derive_var_extreme_dtm()], [derive_var_extreme_dt()]
#'
#' @family deprecated
#' @keywords deprecated
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
#'   date = convert_dtc_to_dt(LBDTC)
#' )
#'
#' # death date from ADSL including traceability variables
#' death_date <- date_source(
#'   dataset_name = "adsl",
#'   date = DTHDT,
#'   set_values_to = exprs(
#'     LALVDOM = "ADSL",
#'     LALVVAR = "DTHDT"
#'   )
#' )
date_source <- function(dataset_name,
                        filter = NULL,
                        date,
                        set_values_to = NULL) {
  deprecate_inform(
    when = "1.2.0",
    what = "date_source()",
    with = "event()",
    details = c(
      x = "This message will turn into a warning at the beginning of 2026.",
      i = "See admiral's deprecation guidance:
      https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )

  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enexpr(filter), optional = TRUE),
    date = assert_expr(enexpr(date)),
    set_values_to = assert_expr_list(set_values_to, named = TRUE, optional = TRUE)
  )
  class(out) <- c("date_source", "source", "list")
  out
}
