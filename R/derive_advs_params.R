#' Adds a Parameter for Mean Arterial Pressure
#'
#' Adds a record for mean arterial pressure (MAP) for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{2DIABP + SYSBP}{3}}{(2DIABP + SYSBP) / 3}
#' if it is based on diastolic and systolic blood pressure and
#' \deqn{DIABP + 0.01 e^{4.14 - \frac{40.74}{HR}} (SYSBP - DIABP)}{
#' DIABP + 0.01 exp(4.14 - 40.74 / HR) (SYSBP - DIABP)}
#' if it is based on diastolic, systolic blood pressure, and heart rate.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and the `unit_var` parameter,
#'   `PARAMCD`, and `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `sysbp_code`, `diabp_code`
#'   and `hr_code`.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param sysbp_code Systolic blood pressure parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the systolic blood pressure assessments.
#'
#'   Permitted Values: character value
#'
#' @param diabp_code Diastolic blood pressure parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the diastolic blood pressure assessments.
#'
#'   Permitted Values: character value
#'
#' @param hr_code Heart rate parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the heart rate assessments.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param unit_var Variable providing the unit of the parameter
#'
#'   For the new parameter the variable is set to the value of the variable for
#'   systolic blood pressure.
#'
#'   Permitted Values: A variable of the input dataset
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation advs
#'
#' @export
#'
#' @examples
#' # Derive MAP based on diastolic and systolic blood pressure
#' advs <- tibble::tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "mmHg", "BASELINE",
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "mmHg", "WEEK 2",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "BASELINE",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "WEEK 2",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "mmHg", "BASELINE",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "mmHg", "WEEK 2",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "mmHg", "BASELINE",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "mmHg", "WEEK 2"
#' )
#'
#' derive_param_map(advs,
#'                  by_vars = vars(USUBJID, VISIT),
#'                  unit_var = AVALU,
#'                  set_values_to = vars(PARAM = "Mean Arterial Pressure (mmHg)"))
#'
#' # Derive MAP based on diastolic and systolic blood pressure and heart rate
#' advs <- tibble::tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU,      ~VISIT,
#'   "01-701-1015", "PULSE",  "Pulse (beats/min)"              ,  59,   "beats/min", "BASELINE",
#'   "01-701-1015", "PULSE",  "Pulse (beats/min)"              ,  61,   "beats/min", "WEEK 2",
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "mmHg",      "BASELINE",
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "mmHg",      "WEEK 2",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg",      "BASELINE",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg",      "WEEK 2",
#'   "01-701-1028", "PULSE",  "Pulse (beats/min)"              ,  62,   "beats/min", "BASELINE",
#'   "01-701-1028", "PULSE",  "Pulse (beats/min)"              ,  77,   "beats/min", "WEEK 2",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "mmHg",      "BASELINE",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "mmHg",      "WEEK 2",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "mmHg",      "BASELINE",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "mmHg",      "WEEK 2"
#' )
#'
#' derive_param_map(advs,
#'   hr_code = "PULSE",
#'   by_vars = vars(USUBJID, VISIT),
#'   unit_var = AVALU,
#'   set_values_to = vars(PARAM = "Mean Arterial Pressure (mmHg)"))

derive_param_map <- function(dataset,
                             filter = NULL,
                             new_param = "MAP",
                             sysbp_code = "SYSBP",
                             diabp_code = "DIABP",
                             hr_code = NULL,
                             by_vars,
                             unit_var = NULL,
                             set_values_to = NULL,
                             drop_values_from = NULL) {
  assert_character_scalar(new_param)
  assert_character_scalar(sysbp_code)
  assert_character_scalar(diabp_code)
  assert_character_scalar(hr_code, optional = TRUE)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset,
                    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var))
  assert_param_does_not_exist(dataset, new_param)
  assert_varval_list(set_values_to, optional = TRUE)

  if (!quo_is_null(unit_var)) {
    unit <-
      unique(filter(dataset, PARAMCD == sysbp_code &
                      !is.na(!!unit_var))[[as_string(quo_get_expr(unit_var))]])
    set_unit_var <- vars(!!unit_var := unit)
  }
  else {
    set_unit_var <- NULL
  }
  if (is.null(hr_code)) {
    analysis_value <-
      expr((2 * !!sym(paste0("AVAL.", diabp_code)) +!!sym(paste0("AVAL.", sysbp_code))) / 3)
  }
  else {
    analysis_value <-
      expr(!!sym(paste0("AVAL.", diabp_code)) +
             0.01 * exp(4.14 - 40.74 / !!sym(paste0("AVAL.", hr_code))) *
             (!!sym(paste0("AVAL.", sysbp_code)) - !!sym(paste0("AVAL.", diabp_code))))
  }
  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(sysbp_code, diabp_code, hr_code),
    by_vars = by_vars,
    analysis_value = !!analysis_value,
    set_values_to = vars(PARAMCD = !!new_param,
                         !!!set_unit_var,
                         !!!set_values_to),
    drop_values_from = drop_values_from
  )
}
