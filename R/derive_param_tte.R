#' Derive a Time-to-event Parameter
#'
#' Add a time-to-event parameter to the input dataset.
#'
#' @param dataset Input dataset
#'
#'   The `USUBJID` variable is expected.
#'
#' @param dataset_adsl ADSL input dataset
#'   The variables specified for `start_date` and `start_imputation_flag` are
#'   expected.
#'
#' @param start_date Time to event origin date
#'   The variable `STARTDT` is set to the specified date. The value is taken
#'   from the ADSL dataset.
#'
#'   If the event or censoring date is before the origin date, `ADT` is set to
#'   the origin date.
#'
#'   If the specified variable is imputed, the corresponding date imputation
#'   flag must specified for `start_imputation_flag`.
#'
#' @param start_imputation_flag
#'   If the start date is imputed, the corresponding date imputation flag must
#'   be specified. The variable `STARTDTF` is set to the specified variable.
#'
#' @param event_conditions Sources and conditions defining events. A list of
#'   `tte_source()` objects is expected.
#'
#' @param censor_conditions Sources and conditions defining censorings. A list
#'   of `tte_source()` objects is expected.
#'
#' @param set_values_to A named list returned by `vars()` defining the variables
#'   to be set for the new parameter, e.g. `vars(PARAMCD = "OS", PARAM =
#'   "Overall Survival")`. The values must be a symbol, a character string, a
#'   numeric value, or `NA`.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
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
#'   the datepart is copied. If the source variable is a character variable, it
#'   is converted to a date. If the date is incomplete, it is imputed as the
#'   first possible date.
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
#'   the datepart is copied. If the source variable is a character variable, it
#'   is converted to a date. If the date is incomplete, it is imputed as the
#'   first possible date.
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
#'   observation is added to the output dataset. If an event is available the
#'   event observations is added. Otherwise the censoring observation is added.
#'
#'   Finally the variables as defined by the `set_values_to` parameter are added
#'   for the new observations.
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation adtte
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("adsl")
#'
#' death <- tte_source(
#'   dataset = adsl,
#'   filter = DTHFL == "Y",
#'   date = DTHDT,
#'   set_values_to =vars(
#'     EVENTDESC = "DEATH",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "DTHDT"))
#'
#' lstalv <- tte_source(
#'   dataset = adsl,
#'   date = LSTALVDT,
#'   censor = 1,
#'   set_values_to = vars(
#'     EVENTDESC = "LAST KNOWN ALIVE DATE",
#'     SRCDOM = "ADSL",
#'     SRCVAR = "LSTALVDT"))
#'
#' derive_param_tte(
#'   event_conditions = list(death),
#'   censor_conditions = list(lstalv),
#'   set_values_to = vars(PARAMCD = "OS",
#'                        PARAM = "Overall Survival"))
derive_param_tte <- function(dataset = NULL,
                             dataset_adsl,
                             start_date = TRTSDT,
                             start_imputation_flag = NULL,
                             event_conditions,
                             censor_conditions,
                             set_values_to,
                             subject_keys = vars(STUDYID, USUBJID)) {
  assert_data_frame(dataset, optional = TRUE)
  start_date <- assert_symbol(enquo(start_date))
  start_imputation_flag <- assert_symbol(enquo(start_imputation_flag),
                                         optional = TRUE)
  assert_data_frame(dataset_adsl,
                    required_vars = quo_c(start_date, start_imputation_flag))
  assert_vars(subject_keys)
  assert_list_of(event_conditions, "tte_source")
  assert_list_of(censor_conditions, "tte_source")
  assert_varval_list(set_values_to, optional = TRUE)

  event_data <- filter_date_sources(sources = event_conditions,
                                    subject_keys = subject_keys,
                                    mode = "first") %>%
    mutate(temp_event = 1)
  censor_data <- filter_date_sources(sources = censor_conditions,
                                     subject_keys = subject_keys,
                                     mode = "last") %>%
    mutate(temp_event = 0)

  adsl_vars = vars(!!!subject_keys,
                   STARTDT = !!start_date)
  if (!quo_is_null(start_imputation_flag)) {
    adls_vars = vars(!!!adsl_vars,
                     STARTDTF = !!start_imputation_flag)
  }
  adsl <- dataset_adsl %>%
    select(!!!adsl_vars)

  new_param <- filter_extreme(
    bind_rows(event_data, censor_data),
    by_vars = subject_keys,
    order = vars(temp_event),
    mode = "last") %>%
    left_join(adsl,
              by = vars2chr(subject_keys)) %>%
    mutate(!!!set_values_to,
           ADT = max(ADT, STARTDT)) %>%
    select(-starts_with("temp_"))

  bind_rows(dataset, new_param)
}

filter_date_sources <- function(sources,
                                subject_keys,
                                mode) {
  assert_vars(subject_keys)

  assert_list_of(sources, "tte_source")

  data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    date <- quo_get_expr(sources[[i]]$date)
    data[[i]] <- sources[[i]]$dataset %>%
      filter_if(sources[[i]]$filter) %>%
      filter_extreme(
        order = vars(!!date),
        by_vars = subject_keys,
        mode = mode,
        check_type = "none"
      )
    # add date variable and accompanying variables
    if (is.Date(data[[i]][[as_string(date)]])) {
      data[[i]] <- transmute(data[[i]],
                             !!!subject_keys,
                             !!!sources[[i]]$set_values_to,
                             CNSR = sources[[i]]$censor,
                             ADT = !!date)
    } else if (is.instant(data[[i]][[as_string(date)]])) {
      data[[i]] <- transmute(data[[i]],
                             !!!subject_keys,
                             !!!sources[[i]]$set_values_to,
                             CNSR = sources[[i]]$censor,
                             ADT = date(!!date))
    } else {
      data[[i]] <- transmute(
        data[[i]],
        !!!subject_keys,
        !!!sources[[i]]$set_values_to,
        CNSR = sources[[i]]$censor,
        ADT = convert_dtc_to_dt(!!date,
                                date_imputation = "first")
      )
    }
  }

  # put all source data into one dataset and select first or last date per subject
  data %>%
    bind_rows() %>%
    filter(!is.na(ADT)) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(ADT),
      mode = mode,
      check_type = "none"
    )
}

#' Create an `tte_source` object
#'
#' The `tte_source` object is used to define events and possible censorings as
#' input for the `derive_param_tte()` function.
#'
#' @param dataset The source dataset
#'
#' @param filter An unquoted condition for selecting the observations from
#'   `dataset` which are events or possible censoring time points.
#'
#' @param date A variable providing the date of the event or censoring. A date,
#'   a datetime, or a character variable containing ISO 8601 dates can be
#'   specified. An unquoted symbol is expected.
#'
#' @param censor Censoring value
#'   CDISC strongly recommends using `0` for events and positive integers for
#'   censoring.
#'
#' @param set_values_to A named list returned by `vars()` defining the
#'   variables to be set for the event or censoring, e.g. `vars(EVENTDESC =
#'   “DEATH”, SRCDOM = “ADSL”, SRCVAR = “DTHDT”)`. The values must be a symbol,
#'   a character string, a numerig value, or `NA`.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "tte_source".
tte_source <- function(dataset,
                       filter = NULL,
                       date,
                       censor = 0,
                       set_values_to = NULL) {
  out <- list(
    dataset = assert_data_frame(dataset),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    censor = assert_integer_scalar(censor),
    set_values_to = assert_varval_list(set_values_to, optional = TRUE)
  )
  class(out) <- c("tte_source", "list")
  out
}
