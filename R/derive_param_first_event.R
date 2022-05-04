#' Add a First Event Parameter
#'
#' Add a new parameter for the first event occurring within a parameter of the
#' input dataset. `AVALC` and `AVAL` indicate if an event occurred and `ADT` is
#' set to the date of the first event. For example, the function can derive a
#' parameter for the first disease progression.
#'
#' @param dataset Input dataset
#'
#'   The variables `PARAMCD`, `ADT`, and those specified by the `subject_keys`
#'   parameter are expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `subject_keys` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param source_param Source parameter
#'
#'   Only the observations from the input dataset where `PARAMCD` equals the
#'   specified value are considered for looking for an event.
#'
#' @param condition Event condition
#'
#'   For subjects with at least one observation where the condition is fulfilled
#'   `AVALC` is set to `"Y"`, `AVAL` to `1`, and `ADT` to the first date where
#'   the condition is fulfilled.
#'
#'   For all other subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to
#'   `NA`.
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
#'   1. The input dataset is restricted to observations where `PARAMCD` equals
#'   `source_param`.
#'   1. For each subject (with respect to the variables specified for the
#'   `subject_keys` parameter) the first observation (with respect `ADT`) where
#'   the event condition (`condition` parameter) is fulfilled is selected.
#'   1. For each observation in `dataset_adsl` a new observation is created. For
#'   subjects with event `AVALC` is set to `"Y"`, `AVAL` to `1`, and `ADT` to
#'   the first date where the event condition is fulfilled. For all other
#'   subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to `NA`.
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
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3"
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
#'   source_param = "OVR",
#'   condition = AVALC == "PD",
#'   set_values_to = vars(
#'     PARAMCD = "PD",
#'     PARAM = "Disease Progression",
#'     ANL01FL = "Y"
#'   )
#' )
derive_param_first_event <- function(dataset,
                                     dataset_adsl,
                                     source_param,
                                     condition,
                                     subject_keys = vars(STUDYID, USUBJID),
                                     set_values_to,
                                     check_type = "warning") {
  # Check input parameters
  assert_character_scalar(source_param)
  condition <- assert_filter_cond(enquo(condition))
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = vars(!!!subject_keys, PARAMCD, ADT))
  assert_data_frame(dataset_adsl, required_vars = subject_keys)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_varval_list(set_values_to)

  # Create new observations
  new_obs <- derive_vars_merged(
    dataset_adsl,
    dataset_add = filter(dataset, PARAMCD == source_param),
    filter_add = !!condition,
    by_vars = subject_keys,
    order = vars(ADT),
    new_vars = vars(ADT),
    mode = "first",
    check_type = check_type
  ) %>%
    mutate(
      AVALC = if_else(!is.na(ADT), "Y", "N"),
      AVALN = if_else(!is.na(ADT), 1, 0),
      !!!set_values_to
    )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
