#' Add a First Event Parameter
#'
#' Add a new parameter for the first event occurring in a dataset. `AVALC` and
#' `AVAL` indicate if an event occurred and `ADT` is set to the date of the
#' first event. For example, the function can derive a parameter for the first
#' disease progression.
#'
#' @param dataset Input dataset
#'
#'   The `PARAMCD` variable is expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `subject_keys` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param dataset_source Source dataset
#'
#'   All observations in the specified dataset fulfilling the condition
#'   specified by `filter_source` are considered as event.
#'
#'   The variables specified by the `subject_keys` and
#'   `date_var` parameter are expected.
#'
#' @param filter_source Source filter
#'
#'   All observations in `dataset_source` fulfilling the specified condition are
#'   considered as event.
#'
#'   For subjects with at least one event `AVALC` is set to `"Y"`, `AVAL` to
#'   `1`, and `ADT` to the first date where the condition is fulfilled.
#'
#'   For all other subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to
#'   `NA`.
#'
#' @param date_var Date variable
#'
#'   Date variable in the source dataset (`dataset_source`). The variable is
#'   used to sort the source dataset. `ADT` is set to the specified variable for
#'   events.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `vars()` defining the variables to be set for the
#'   new parameter, e.g. `vars(PARAMCD = "PD", PARAM = "Disease Progression")`
#'   is expected. The values must be symbols, character strings, numeric values,
#'   or `NA`.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `vars()` is expected.
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the
#'   observations of the input dataset restricted to the source parameter
#'   (`source_param`) are not unique with respect to the subject keys
#'   (`subject_key` parameter) and `ADT`.
#'
#'   *Default*: `"warning"`
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details
#'   1. The input dataset is restricted to observations fulfilling
#'   `filter_source`.
#'   1. For each subject (with respect to the variables specified for the
#'   `subject_keys` parameter) the first observation (with respect to
#'   `date_var`) where the event condition (`filter_source` parameter) is
#'   fulfilled is selected.
#'   1. For each observation in `dataset_adsl` a new observation is created. For
#'   subjects with event `AVALC` is set to `"Y"`, `AVAL` to `1`, and `ADT` to
#'   the first date where the event condition is fulfilled. For all other
#'   subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to `NA`.
#'   For subjects with event all variables from `dataset_source` are kept. For
#'   subjects without event all variables which are in both `dataset_adsl` and
#'   `dataset_source` are kept.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the new observations.
#'   1. The new observations are added to input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with a new parameter indicating if and when an
#'   event occurred
#'
#' @keywords derivation bds
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(lubridate)
#'
#' # Derive a new parameter for the first disease progression (PD)
#' adsl <- tibble::tribble(
#'   ~USUBJID, ~DTHDT,
#'   "1",      ymd("2022-05-13"),
#'   "2",      ymd(""),
#'   "3",      ymd("")
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' adrs <- tibble::tribble(
#'   ~USUBJID, ~ADTC,        ~AVALC,
#'   "1",      "2020-01-02", "PR",
#'   "1",      "2020-02-01", "CR",
#'   "1",      "2020-03-01", "CR",
#'   "1",      "2020-04-01", "SD",
#'   "2",      "2021-06-15", "SD",
#'   "2",      "2021-07-16", "PD",
#'   "2",      "2021-09-14", "PD"
#' ) %>%
#'   mutate(
#'     STUDYID = "XX1234",
#'     ADT = ymd(ADTC),
#'     PARAMCD = "OVR",
#'     PARAM = "Overall Response",
#'     ANL01FL = "Y"
#'   ) %>%
#'   select(-ADTC)
#'
#' derive_param_first_event(
#'   adrs,
#'   dataset_adsl = adsl,
#'   dataset_source = adrs,
#'   filter_source = PARAMCD == "OVR" & AVALC == "PD",
#'   date_var = ADT,
#'   set_values_to = vars(
#'     PARAMCD = "PD",
#'     PARAM = "Disease Progression",
#'     ANL01FL = "Y"
#'   )
#' )
#'
#' # derive parameter indicating death
#' derive_param_first_event(
#'   dataset = adrs,
#'   dataset_adsl = adsl,
#'   dataset_source = adsl,
#'   filter_source = !is.na(DTHDT),
#'   date_var = DTHDT,
#'   set_values_to = vars(
#'     PARAMCD = "DEATH",
#'     PARAM = "Death",
#'     ANL01FL = "Y"
#'   )
#' )
derive_param_first_event <- function(dataset,
                                     dataset_adsl,
                                     dataset_source,
                                     filter_source,
                                     date_var,
                                     subject_keys = vars(STUDYID, USUBJID),
                                     set_values_to,
                                     check_type = "warning") {
  # Check input parameters
  filter_source <- assert_filter_cond(enquo(filter_source))
  date_var <- assert_symbol(enquo(date_var))
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = vars(PARAMCD))
  assert_data_frame(dataset_source, required_vars = vars(!!!subject_keys, !!date_var))
  assert_data_frame(dataset_adsl, required_vars = subject_keys)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  # Create new observations
  source_vars <- colnames(dataset_source)
  adsl_vars <- colnames(dataset_adsl)

  events <- dataset_source %>%
    filter_if(filter_source) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(!!date_var),
      mode = "first",
      check_type = check_type
  )
  noevents <- anti_join(
    select(dataset_adsl, intersect(source_vars, adsl_vars)),
    select(events, !!!subject_keys))
  new_obs <- bind_rows(events, noevents) %>%
    mutate(
      ADT = !!date_var,
      AVALC = if_else(!is.na(ADT), "Y", "N"),
      AVAL = if_else(!is.na(ADT), 1, 0),
      !!!set_values_to
    )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
