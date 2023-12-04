#' Derive a Time-to-Event Parameter
#'
#' Add a time-to-event parameter to the input dataset.
#'
#' @param dataset `r roxygen_param_dataset()`
#'   `PARAMCD` is expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `start_date`, and
#'   `subject_keys` are expected.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `tte_source()` refers to the dataset provided in the list.
#'
#' @param by_vars By variables
#'
#'   If the parameter is specified, separate time to event parameters are
#'   derived for each by group.
#'
#'   The by variables must be in at least one of the source datasets. Each
#'   source dataset must contain either all by variables or none of the by
#'   variables.
#'
#'   The by variables are not included in the output dataset.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param start_date Time to event origin date
#'
#'   The variable `STARTDT` is set to the specified date. The value is taken
#'   from the ADSL dataset.
#'
#'   If the event or censoring date is before the origin date, `ADT` is set to
#'   the origin date.
#'
#' @param event_conditions Sources and conditions defining events
#'
#'   A list of `event_source()` objects is expected.
#'
#' @param censor_conditions Sources and conditions defining censorings
#'
#'   A list of `censor_source()` objects is expected.
#'
#' @param create_datetime Create datetime variables?
#'
#'   If set to `TRUE`, variables `ADTM` and `STARTDTM` are created. Otherwise,
#'   variables `ADT` and `STARTDT` are created.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `exprs()` defining the variables to be set for the
#'   new parameter, e.g. `exprs(PARAMCD = "OS", PARAM = "Overall Survival")` is
#'   expected. The values must be symbols, character strings, numeric values,
#'   expressions, or `NA`.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `exprs()` is expected.
#'
#' @details The following steps are performed to create the observations of the
#'   new parameter:
#'
#'   **Deriving the events:**
#'
#'   \enumerate{ \item For each event source dataset the observations as
#'   specified by the `filter` element are selected. Then for each patient the
#'   first observation (with respect to `date`) is selected.
#'
#'   \item The `ADT` variable is set to the variable specified by the
#'   \code{date} element. If the date variable is a datetime variable, only
#'   the datepart is copied.
#'
#'   \item The `CNSR` variable is added and set to the \code{censor} element.
#'
#'   \item The variables specified by the \code{set_values_to} element are
#'   added.
#'
#'   \item The selected observations of all event source datasets are combined into a
#'   single dataset.
#'
#'   \item For each patient the first observation (with respect to the `ADT`
#'   variable) from the single dataset is selected. }
#'
#'   **Deriving the censoring observations:**
#'
#'   \enumerate{ \item For each censoring source dataset the observations as
#'   specified by the `filter` element are selected. Then for each patient the
#'   last observation (with respect to `date`) is selected.
#'
#'   \item The `ADT` variable is set to the variable specified by the
#'   \code{date} element. If the date variable is a datetime variable, only
#'   the datepart is copied.
#'
#'   \item The `CNSR` variable is added and set to the \code{censor} element.
#'
#'   \item The variables specified by the \code{set_values_to} element are
#'   added.
#'
#'   \item The selected observations of all censoring source datasets are
#'   combined into a single dataset.
#'
#'   \item For each patient the last observation (with respect to the `ADT`
#'   variable) from the single dataset is selected. }
#'
#'   For each subject (as defined by the `subject_keys` parameter) an
#'   observation is selected. If an event is available, the event observation is
#'   selected. Otherwise the censoring observation is selected.
#'
#'   Finally:
#'   1. The variable specified for `start_date` is joined from the
#'   ADSL dataset. Only subjects in both datasets are kept,
#'   i.e., subjects with both an event or censoring and an observation in
#'   `dataset_adsl`.
#'   1. The variables as defined by the `set_values_to` parameter are added.
#'   1. The `ADT`/`ADTM` variable is set to the maximum of `ADT`/`ADTM` and
#'   `STARTDT`/`STARTDTM` (depending on the `create_datetime` parameter).
#'   1. The new observations are added to the output dataset.
#'
#'
#' @return The input dataset with the new parameter added
#'
#' @family der_prm_tte
#' @keywords der_prm_tte
#'
#' @seealso [event_source()], [censor_source()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#' data("admiral_adsl")
#'
#' adsl <- admiral_adsl
#'
#' # derive overall survival parameter
#' death <- event_source(
#'   dataset_name = "adsl",
#'   filter = DTHFL == "Y",
#'   date = DTHDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "DEATH",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "DTHDT"
#'   )
#' )
#'
#' last_alive_dt <- censor_source(
#'   dataset_name = "adsl",
#'   date = LSTALVDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "LAST DATE KNOWN ALIVE",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "LSTALVDT"
#'   )
#' )
#'
#' derive_param_tte(
#'   dataset_adsl = adsl,
#'   event_conditions = list(death),
#'   censor_conditions = list(last_alive_dt),
#'   source_datasets = list(adsl = adsl),
#'   set_values_to = exprs(
#'     PARAMCD = "OS",
#'     PARAM = "Overall Survival"
#'   )
#' ) %>%
#'   select(-STUDYID) %>%
#'   filter(row_number() %in% 20:30)
#'
#' # derive duration of response
#' # only observations for subjects in dataset_adsl are created
#' adsl <- tribble(
#'   ~USUBJID, ~DTHFL, ~DTHDT,            ~RSPDT,
#'   "01",     "Y",    ymd("2021-06-12"), ymd("2021-03-04"),
#'   "02",     "N",    NA,                NA,
#'   "03",     "Y",    ymd("2021-08-21"), NA,
#'   "04",     "N",    NA,                ymd("2021-04-14")
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' adrs <- tribble(
#'   ~USUBJID, ~AVALC, ~ADT,              ~ASEQ,
#'   "01",     "SD",   ymd("2021-01-03"), 1,
#'   "01",     "PR",   ymd("2021-03-04"), 2,
#'   "01",     "PD",   ymd("2021-05-05"), 3,
#'   "02",     "PD",   ymd("2021-02-03"), 1,
#'   "04",     "SD",   ymd("2021-02-13"), 1,
#'   "04",     "PR",   ymd("2021-04-14"), 2,
#'   "04",     "CR",   ymd("2021-05-15"), 3
#' ) %>%
#'   mutate(STUDYID = "AB42", PARAMCD = "OVR")
#'
#' pd <- event_source(
#'   dataset_name = "adrs",
#'   filter = AVALC == "PD",
#'   date = ADT,
#'   set_values_to = exprs(
#'     EVENTDESC = "PD",
#'     SRCDOM = "ADRS",
#'     SRCVAR = "ADTM",
#'     SRCSEQ = ASEQ
#'   )
#' )
#'
#' death <- event_source(
#'   dataset_name = "adsl",
#'   filter = DTHFL == "Y",
#'   date = DTHDT,
#'   set_values_to = exprs(
#'     EVENTDESC = "DEATH",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "DTHDT"
#'   )
#' )
#'
#' lastvisit <- censor_source(
#'   dataset_name = "adrs",
#'   date = ADT,
#'   censor = 1,
#'   set_values_to = exprs(
#'     EVENTDESC = "LAST TUMOR ASSESSMENT",
#'     SRCDOM = "ADRS",
#'     SRCVAR = "ADTM",
#'     SRCSEQ = ASEQ
#'   )
#' )
#'
#' first_response <- censor_source(
#'   dataset_name = "adsl",
#'   date = RSPDT,
#'   censor = 1,
#'   set_values_to = exprs(
#'     EVENTDESC = "FIRST RESPONSE",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "RSPDT"
#'   )
#' )
#'
#' derive_param_tte(
#'   dataset_adsl = filter(adsl, !is.na(RSPDT)),
#'   start_date = RSPDT,
#'   event_conditions = list(pd, death),
#'   censor_conditions = list(lastvisit, first_response),
#'   source_datasets = list(adsl = adsl, adrs = adrs),
#'   set_values_to = exprs(
#'     PARAMCD = "DURRSP",
#'     PARAM = "Duration of Response"
#'   )
#' )
#'
#' # derive time to adverse event for each preferred term
#' adsl <- tribble(
#'   ~USUBJID, ~TRTSDT,           ~EOSDT,
#'   "01",     ymd("2020-12-06"), ymd("2021-03-06"),
#'   "02",     ymd("2021-01-16"), ymd("2021-02-03")
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' ae <- tribble(
#'   ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
#'   "01",     "2021-01-03T10:56", 1,      "Flu",
#'   "01",     "2021-03-04",       2,      "Cough",
#'   "01",     "2021",             3,      "Flu"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' ae_ext <- derive_vars_dt(
#'   ae,
#'   dtc = AESTDTC,
#'   new_vars_prefix = "AEST",
#'   highest_imputation = "M",
#'   flag_imputation = "none"
#' )
#'
#' ttae <- event_source(
#'   dataset_name = "ae",
#'   date = AESTDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "AE",
#'     SRCDOM = "AE",
#'     SRCVAR = "AESTDTC",
#'     SRCSEQ = AESEQ
#'   )
#' )
#'
#' eos <- censor_source(
#'   dataset_name = "adsl",
#'   date = EOSDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "END OF STUDY",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "EOSDT"
#'   )
#' )
#'
#' derive_param_tte(
#'   dataset_adsl = adsl,
#'   by_vars = exprs(AEDECOD),
#'   start_date = TRTSDT,
#'   event_conditions = list(ttae),
#'   censor_conditions = list(eos),
#'   source_datasets = list(adsl = adsl, ae = ae_ext),
#'   set_values_to = exprs(
#'     PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
#'     PARAM = paste("Time to First", AEDECOD, "Adverse Event"),
#'     PARCAT1 = "TTAE",
#'     PARCAT2 = AEDECOD
#'   )
#' ) %>%
#'   select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
derive_param_tte <- function(dataset = NULL,
                             dataset_adsl,
                             source_datasets,
                             by_vars = NULL,
                             start_date = TRTSDT,
                             event_conditions,
                             censor_conditions,
                             create_datetime = FALSE,
                             set_values_to,
                             subject_keys = get_admiral_option("subject_keys")) {
  # checking and quoting #
  assert_data_frame(dataset, optional = TRUE)
  assert_vars(by_vars, optional = TRUE)
  start_date <- assert_symbol(enexpr(start_date))
  assert_data_frame(dataset_adsl, required_vars = exprs(!!start_date))
  assert_vars(subject_keys)
  assert_list_of(event_conditions, "event_source")
  assert_list_of(censor_conditions, "censor_source")
  assert_list_of(source_datasets, "data.frame")
  source_names <- names(source_datasets)
  assert_list_element(
    list = event_conditions,
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
  assert_list_element(
    list = censor_conditions,
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
  assert_logical_scalar(create_datetime)
  assert_varval_list(set_values_to, optional = TRUE)
  if (!is.null(set_values_to$PARAMCD) && !is.null(dataset)) {
    assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  }
  if (!is.null(by_vars)) {
    source_datasets <- extend_source_datasets(
      source_datasets = source_datasets,
      by_vars = by_vars
    )
  }

  tmp_event <- get_new_tmp_var(dataset)
  # determine events #
  event_data <- filter_date_sources(
    sources = event_conditions,
    source_datasets = source_datasets,
    by_vars = by_vars,
    create_datetime = create_datetime,
    subject_keys = subject_keys,
    mode = "first"
  ) %>%
    mutate(!!tmp_event := 1L)

  # determine censoring observations #
  censor_data <- filter_date_sources(
    sources = censor_conditions,
    source_datasets = source_datasets,
    by_vars = by_vars,
    create_datetime = create_datetime,
    subject_keys = subject_keys,
    mode = "last"
  ) %>%
    mutate(!!tmp_event := 0L)

  # determine variable to add from ADSL #
  if (create_datetime) {
    date_var <- sym("ADTM")
    start_var <- sym("STARTDTM")
  } else {
    date_var <- sym("ADT")
    start_var <- sym("STARTDT")
  }
  adsl_vars <- exprs(
    !!!subject_keys,
    !!start_var := !!start_date
  )

  start_date_imputation_flag <- gsub("(DT|DTM)$", "DTF", as_name(start_date))
  if (start_date_imputation_flag %in% colnames(dataset_adsl) &&
    as_name(start_date) != start_date_imputation_flag) {
    adsl_vars <- exprs(
      !!!adsl_vars,
      STARTDTF = !!sym(start_date_imputation_flag)
    )
  }

  start_time_imputation_flag <- gsub("DTM$", "TMF", as_name(start_date))
  if (start_time_imputation_flag %in% colnames(dataset_adsl) &&
    as_name(start_date) != start_time_imputation_flag) {
    adsl_vars <- exprs(
      !!!adsl_vars,
      STARTTMF = !!sym(start_time_imputation_flag)
    )
  }

  adsl <- dataset_adsl %>%
    select(!!!adsl_vars)

  # create observations for new parameter #
  new_param <- filter_extreme(
    bind_rows(event_data, censor_data),
    by_vars = expr_c(subject_keys, by_vars),
    order = exprs(!!tmp_event),
    mode = "last"
  ) %>%
    inner_join(
      adsl,
      by = vars2chr(subject_keys)
    ) %>%
    process_set_values_to(set_values_to) %>%
    mutate(!!date_var := pmax(!!date_var, !!start_var, na.rm = TRUE)) %>%
    remove_tmp_vars()

  if (!is.null(by_vars)) {
    if (!is.null(set_values_to$PARAMCD)) {
      assert_one_to_one(new_param, exprs(PARAMCD), by_vars)
    }

    # -vars2chr(by_vars) does not work for 3.5 #
    new_param <- select(new_param, !!!negate_vars(by_vars))
  }

  # add new parameter to input dataset #
  bind_rows(dataset, new_param)
}

#' Select the First or Last Date from Several Sources
#'
#' Select for each subject the first or last observation with respect to a date
#' from a list of sources.
#'
#' @param sources Sources
#'
#'    A list of `tte_source()` objects is expected.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `tte_source()` refers to the dataset provided in the list.
#'
#' @param by_vars By variables
#'
#'   If the parameter is specified, for each by group the observations are
#'   selected separately.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param create_datetime Create datetime variable?
#'
#'   If set to `TRUE`, variables `ADTM` is created. Otherwise, variables `ADT`
#'   is created.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `exprs()` is expected.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each subject the first observation with
#'   respect to the date is included in the output dataset. If `"last"` is
#'   specified, the last observation is included in the output dataset.
#'
#'   Permitted Values:  `"first"`, `"last"`
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{ \item For each source dataset the observations as specified by
#'   the `filter` element are selected. Then for each patient the first or last
#'   observation (with respect to `date`) is selected.
#'
#'   \item The `ADT` variable is set to the variable specified by the
#'   \code{date} element. If the date variable is a datetime variable, only
#'   the datepart is copied. If the source variable is a character variable, it
#'   is converted to a date. If the date is incomplete, it is imputed as
#'   the first possible date.
#'
#'   \item The `CNSR` is added and set to the value of the \code{censor}
#'   element.
#'
#'   \item The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   \item For each patient the first or last observation (with respect to the
#'   `ADT` variable) from the single dataset is selected. }
#'
#' @return A dataset with one observation per subject as described in the
#'   "Details" section.
#'
#' @noRd
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adsl <- tribble(
#'   ~USUBJID, ~TRTSDT,           ~EOSDT,
#'   "01",     ymd("2020-12-06"), ymd("2021-03-06"),
#'   "02",     ymd("2021-01-16"), ymd("2021-02-03")
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' ae <- tribble(
#'   ~USUBJID, ~AESTDTC,     ~AESEQ, ~AEDECOD,
#'   "01",     "2021-01-03", 1,      "Flu",
#'   "01",     "2021-03-04", 2,      "Cough",
#'   "01",     "2021-01-01", 3,      "Flu"
#' ) %>%
#'   mutate(
#'     STUDYID = "AB42",
#'     AESTDT = ymd(AESTDTC)
#'   )
#'
#' ttae <- event_source(
#'   dataset_name = "ae",
#'   date = AESTDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "AE",
#'     SRCDOM = "AE",
#'     SRCVAR = "AESTDTC",
#'     SRCSEQ = AESEQ
#'   )
#' )
#'
#' filter_date_sources(
#'   sources = list(ttae),
#'   source_datasets = list(adsl = adsl, ae = ae),
#'   by_vars = exprs(AEDECOD),
#'   create_datetime = FALSE,
#'   subject_keys = get_admiral_option("subject_keys"),
#'   mode = "first"
#' )
filter_date_sources <- function(sources,
                                source_datasets,
                                by_vars,
                                create_datetime = FALSE,
                                subject_keys,
                                mode) {
  assert_list_of(sources, "tte_source")
  assert_list_of(source_datasets, "data.frame")
  assert_logical_scalar(create_datetime)
  assert_vars(subject_keys)
  assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE
  )

  if (create_datetime) {
    date_var <- sym("ADTM")
  } else {
    date_var <- sym("ADT")
  }

  data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    source_date <- sources[[i]]$date
    source_dataset <- source_datasets[[sources[[i]]$dataset_name]]
    if (is.symbol(source_date)) {
      source_date_var <- source_date
    } else {
      source_date_var <- get_new_tmp_var(dataset = source_dataset, prefix = "tmp_date")
      source_dataset <- mutate(
        source_dataset,
        !!source_date_var := !!source_date
      )
    }
    assert_date_var(
      dataset = source_dataset,
      var = !!source_date_var,
      dataset_name = sources[[i]]$dataset_name
    )
    data[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      filter_extreme(
        order = exprs(!!source_date_var),
        by_vars = expr_c(subject_keys, by_vars),
        mode = mode,
        check_type = "none"
      )

    # add date variable and accompanying variables
    if (create_datetime) {
      date_derv <- exprs(!!date_var := as_datetime(!!source_date_var))
    } else {
      date_derv <- exprs(!!date_var := date(!!source_date_var))
    }

    data[[i]] <- transmute(
      data[[i]],
      !!!by_vars,
      !!!subject_keys,
      !!!sources[[i]]$set_values_to,
      CNSR = sources[[i]]$censor,
      !!!date_derv
    )
  }

  # put all source data into one dataset and select first or last date per subject
  data %>%
    bind_rows() %>%
    filter(!is.na(!!date_var)) %>%
    filter_extreme(
      by_vars = expr_c(subject_keys, by_vars),
      order = exprs(!!date_var),
      mode = mode,
      check_type = "none"
    )
}

#' Add By Groups to All Datasets if Necessary
#'
#' The function ensures that the by variables are contained in all source
#' datasets.
#'
#' @details
#'   1. The by groups are determined as the union of the by groups occurring in
#'   the source datasets.
#'   1. For all source datasets which do not contain the by variables the source
#'   dataset is replaced by the cartesian product of the source dataset and the
#'   by groups.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. Each dataset must contain either all
#'   by variables or none of the by variables.
#'
#' @param by_vars By variables
#'
#' `r roxygen_param_by_vars()`
#'
#'
#' @return The list of extended source datasets
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adsl <- tribble(
#'   ~USUBJID, ~TRTSDT,           ~EOSDT,
#'   "01",     ymd("2020-12-06"), ymd("2021-03-06"),
#'   "02",     ymd("2021-01-16"), ymd("2021-02-03")
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' ae <- tribble(
#'   ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
#'   "01",     "2021-01-03T10:56", 1,      "Flu",
#'   "01",     "2021-03-04",       2,      "Cough",
#'   "01",     "2021",             3,      "Flu"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' extend_source_datasets(
#'   source_datasets = list(adsl = adsl, ae = ae),
#'   by_vars = exprs(AEDECOD)
#' )
#'
#' @noRd
extend_source_datasets <- function(source_datasets,
                                   by_vars) {
  assert_list_of(source_datasets, "data.frame")
  assert_vars(by_vars)
  # ensure that the by variables are contained in all source datasets #
  by_vars_chr <- vars2chr(by_vars)
  # determine which source datasets need to be extended #
  by_groups <- list()
  extend <- vector("list", length(source_datasets))
  for (i in seq_along(source_datasets)) {
    missing_by_vars <- setdiff(by_vars_chr, names(source_datasets[[i]]))
    if (length(missing_by_vars) == 0) {
      # all by variables are included in the source dataset #
      extend[[i]] <- FALSE
      by_groups <-
        append(
          by_groups,
          list(unique(select(
            source_datasets[[i]], !!!by_vars
          )))
        )
    } else if (!setequal(by_vars_chr, missing_by_vars)) {
      # only some of the by variables are included in the source dataset #
      abort(paste(
        "Only",
        paste(setdiff(by_vars_chr, missing_by_vars), collapse = ", "),
        "are included in source dataset",
        names(source_datasets)[[i]],
        ".\n The source dataset must include all or none of the by variables."
      ))
    } else {
      extend[[i]] <- TRUE
    }
  }
  if (length(by_groups) == 0) {
    abort(paste0(
      "The by variables (",
      paste(by_vars_chr, collapse = ", "),
      ") are not contained in any of the source datasets."
    ))
  }
  # extend source datasets #
  by_groups <- unique(bind_rows(by_groups))
  for (i in seq_along(source_datasets)) {
    if (extend[[i]]) {
      source_datasets[[i]] <- crossing(by_groups, source_datasets[[i]])
    }
  }
  source_datasets
}

#' Create a `tte_source` Object
#'
#' The `tte_source` object is used to define events and possible censorings.
#'
#' @param dataset_name The name of the source dataset
#'
#'   The name refers to the dataset provided by the `source_datasets` parameter
#'   of `derive_param_tte()`.
#'
#' @param filter An unquoted condition for selecting the observations from
#'   `dataset` which are events or possible censoring time points.
#'
#' @param date A variable or expression providing the date of the event or
#'   censoring. A date, or a datetime can be specified. An unquoted symbol or
#'   expression is expected.
#'
#'   Refer to `derive_vars_dt()` or `convert_dtc_to_dt()` to impute and derive a
#'   date from a date character vector to a date object.
#'
#' @param censor Censoring value
#'
#'   CDISC strongly recommends using `0` for events and positive integers for
#'   censoring.
#'
#' @param set_values_to A named list returned by `exprs()` defining the variables
#'   to be set for the event or censoring, e.g. `exprs(EVENTDESC = "DEATH",
#'   SRCDOM = "ADSL", SRCVAR = "DTHDT")`. The values must be a symbol, a
#'   character string, a numeric value, an expression, or `NA`.
#'
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_param_tte()], [censor_source()], [event_source()]
#'
#' @return An object of class `tte_source`
tte_source <- function(dataset_name,
                       filter = NULL,
                       date,
                       censor = 0,
                       set_values_to = NULL) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enexpr(filter), optional = TRUE),
    date = assert_expr(enexpr(date)),
    censor = assert_integer_scalar(censor),
    set_values_to = assert_expr_list(
      set_values_to,
      named = TRUE,
      optional = TRUE
    )
  )
  class(out) <- c("tte_source", "source", "list")
  out
}

#' Create an `event_source` Object
#'
#' @description `event_source` objects are used to define events as input for the
#' `derive_param_tte()` function.
#'
#' **Note:** This is a wrapper function for the more generic `tte_source()`.
#'
#' @inheritParams tte_source
#'
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @seealso [derive_param_tte()], [censor_source()]
#'
#' @export
#'
#' @return An object of class `event_source`, inheriting from class `tte_source`
#'
#' @examples
#' # Death event
#'
#' event_source(
#'   dataset_name = "adsl",
#'   filter = DTHFL == "Y",
#'   date = DTHDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "DEATH",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "DTHDT"
#'   )
#' )
event_source <- function(dataset_name,
                         filter = NULL,
                         date,
                         set_values_to = NULL) {
  out <- tte_source(
    dataset_name = assert_character_scalar(dataset_name),
    filter = !!enexpr(filter),
    date = !!assert_expr(enexpr(date)),
    censor = 0,
    set_values_to = set_values_to
  )
  class(out) <- c("event_source", class(out))
  out
}

#' Create a `censor_source` Object
#'
#' @description `censor_source` objects are used to define censorings as input for the
#' `derive_param_tte()` function.
#'
#' **Note:** This is a wrapper function for the more generic `tte_source()`.
#'
#' @inheritParams tte_source
#'
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @seealso [derive_param_tte()], [event_source()]
#'
#' @export
#'
#' @return An object of class `censor_source`, inheriting from class `tte_source`
#'
#' @examples
#' # Last study date known alive censor
#'
#' censor_source(
#'   dataset_name = "adsl",
#'   date = LSTALVDT,
#'   set_values_to = exprs(
#'     EVNTDESC = "ALIVE",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "LSTALVDT"
#'   )
#' )
censor_source <- function(dataset_name,
                          filter = NULL,
                          date,
                          censor = 1,
                          set_values_to = NULL) {
  out <- tte_source(
    dataset_name = assert_character_scalar(dataset_name),
    filter = !!enexpr(filter),
    date = !!assert_expr(enexpr(date)),
    censor = assert_integer_scalar(censor, subset = "positive"),
    set_values_to = set_values_to
  )
  class(out) <- c("censor_source", class(out))
  out
}

#' List all `tte_source` Objects Available in a Package
#'
#' @param package The name of the package in which to search for `tte_source` objects
#'
#' @return
#' A `data.frame` where each row corresponds to one `tte_source` object or `NULL`
#' if `package` does not contain any `tte_source` objects
#'
#'
#' @export
#'
#' @family other_advanced
#' @keywords other_advanced
#'
#' @examples
#' list_tte_source_objects()
list_tte_source_objects <- function(package = "admiral") {
  assert_character_scalar(package)

  if (!requireNamespace(package, quietly = TRUE)) {
    err_msg <- sprintf(
      "No package called '%s' is installed and hence no `tte_source` objects are available",
      package
    )
    abort(err_msg)
  }

  # Get all `tte_source` objects exported by `package`
  exports <- getNamespaceExports(package)
  is_tte_source <- map_lgl(exports, function(obj_name) {
    inherits(getExportedValue(package, obj_name), "tte_source")
  })
  tte_sources <- exports[is_tte_source]

  rows <- lapply(tte_sources, function(obj_name) {
    obj <- getExportedValue(package, obj_name)
    data.frame(
      object = obj_name,
      dataset_name = obj$dataset_name,
      filter = as_label(obj$filter),
      date = as_name(obj$date),
      censor = obj$censor,
      set_values_to = paste(
        paste(
          names(obj$set_values_to),
          map_chr(obj$set_values_to, as_label),
          sep = ": "
        ),
        collapse = "<br>"
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
