#' Adds a Parameter for Mean Arterial Pressure
#'
#' @description Adds a record for mean arterial pressure (MAP) for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' **Note:** This is a wrapper function for the more generic `derive_param_computed()`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `sysbp_code`, `diabp_code`
#'   and `hr_code`.
#'
#' @param sysbp_code Systolic blood pressure parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the systolic blood pressure assessments.
#'
#'   *Permitted Values:* character value
#'
#' @param diabp_code Diastolic blood pressure parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the diastolic blood pressure assessments.
#'
#'   *Permitted Values:* character value
#'
#' @param hr_code Heart rate parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the heart rate assessments.
#'
#'   *Permitted Values:* character value
#'
#' @param set_values_to Variables to be set
#'
#' The specified variables are set to the specified values for the new
#' observations. For example `exprs(PARAMCD = "MAP")` defines the parameter code
#' for the new parameter.
#'
#' *Permitted Values*: List of variable-value pairs
#'
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#' @details
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{2DIABP + SYSBP}{3}}{(2DIABP + SYSBP) / 3}
#' if it is based on diastolic and systolic blood pressure and
#' \deqn{DIABP + 0.01 e^{4.14 - \frac{40.74}{HR}} (SYSBP - DIABP)}{
#' DIABP + 0.01 exp(4.14 - 40.74 / HR) (SYSBP - DIABP)}
#' if it is based on diastolic, systolic blood pressure, and heart rate.
#'
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @seealso [compute_map()]
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
#'   "01-701-1015", "PULSE", "Pulse (beats/min)", 59, "BASELINE",
#'   "01-701-1015", "PULSE", "Pulse (beats/min)", 61, "WEEK 2",
#'   "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
#'   "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "WEEK 2",
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "WEEK 2",
#'   "01-701-1028", "PULSE", "Pulse (beats/min)", 62, "BASELINE",
#'   "01-701-1028", "PULSE", "Pulse (beats/min)", 77, "WEEK 2",
#'   "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "BASELINE",
#'   "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "WEEK 2",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "BASELINE",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "WEEK 2"
#' )
#'
#' # Derive MAP based on diastolic and systolic blood pressure
#' advs %>%
#'   derive_param_map(
#'     by_vars = exprs(USUBJID, VISIT),
#'     set_values_to = exprs(
#'       PARAMCD = "MAP",
#'       PARAM = "Mean Arterial Pressure (mmHg)"
#'     ),
#'     get_unit_expr = extract_unit(PARAM)
#'   ) %>%
#'   filter(PARAMCD != "PULSE")
#'
#' # Derive MAP based on diastolic and systolic blood pressure and heart rate
#' derive_param_map(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   hr_code = "PULSE",
#'   set_values_to = exprs(
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_map <- function(dataset,
                             by_vars,
                             set_values_to = exprs(PARAMCD = "MAP"),
                             sysbp_code = "SYSBP",
                             diabp_code = "DIABP",
                             hr_code = NULL,
                             get_unit_expr,
                             filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, PARAMCD, AVAL))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(sysbp_code)
  assert_character_scalar(diabp_code)
  assert_character_scalar(hr_code, optional = TRUE)
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  assert_unit(dataset, sysbp_code, required_unit = "mmHg", get_unit_expr = !!get_unit_expr)
  assert_unit(dataset, diabp_code, required_unit = "mmHg", get_unit_expr = !!get_unit_expr)

  if (is.null(hr_code)) {
    analysis_value <- expr(
      compute_map(
        diabp = !!sym(paste0("AVAL.", diabp_code)),
        sysbp = !!sym(paste0("AVAL.", sysbp_code))
      )
    )
  } else {
    assert_unit(dataset, hr_code, required_unit = "beats/min", get_unit_expr = !!get_unit_expr)

    analysis_value <- expr(
      compute_map(
        diabp = !!sym(paste0("AVAL.", diabp_code)),
        sysbp = !!sym(paste0("AVAL.", sysbp_code)),
        hr = !!sym(paste0("AVAL.", hr_code))
      )
    )
  }

  withCallingHandlers(
    derive_param_computed(
      dataset,
      filter = !!filter,
      parameters = c(sysbp_code, diabp_code, hr_code),
      by_vars = by_vars,
      set_values_to = exprs(
        AVAL = !!analysis_value,
        !!!set_values_to
      )
    ),
    derive_param_computed_all_na = function(cnd) {
      cli_inform(
        c(
          paste(
            "No computed records were added because for all potential computed",
            "records at least one of the contributing values was {.val {NA}}."
          ),
          "If this is not expected, please check the input data."
        ),
        class = class(cnd)
      )
      cnd_muffle(cnd)
    }
  )
}

#' Compute Mean Arterial Pressure (MAP)
#'
#' Computes mean arterial pressure (MAP) based on diastolic and systolic blood
#' pressure. Optionally heart rate can be used as well.
#'
#' @param diabp Diastolic blood pressure
#'
#'   A numeric vector is expected.
#'
#' @param sysbp Systolic blood pressure
#'
#'   A numeric vector is expected.
#'
#' @param hr Heart rate
#'
#'   A numeric vector or `NULL` is expected.
#'
#'
#' @details
#' \deqn{\frac{2DIABP + SYSBP}{3}}{(2DIABP + SYSBP) / 3}
#' if it is based on diastolic and systolic blood pressure and
#' \deqn{DIABP + 0.01 e^{4.14 - \frac{40.74}{HR}} (SYSBP - DIABP)}{
#' DIABP + 0.01 exp(4.14 - 40.74 / HR) (SYSBP - DIABP)}
#' if it is based on diastolic, systolic blood pressure, and heart rate.
#'
#' Usually this computation function can not be used with `%>%`.
#'
#' @return A numeric vector of MAP values
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @seealso [derive_param_map()]
#'
#' @examples
#' # Compute MAP based on diastolic and systolic blood pressure
#' compute_map(diabp = 51, sysbp = 121)
#'
#' # Compute MAP based on diastolic and systolic blood pressure and heart rate
#' compute_map(diabp = 51, sysbp = 121, hr = 59)
compute_map <- function(diabp, sysbp, hr = NULL) {
  assert_numeric_vector(diabp)
  assert_numeric_vector(sysbp)
  assert_numeric_vector(hr, optional = TRUE)

  if (is.null(hr)) {
    (2 * diabp + sysbp) / 3
  } else {
    diabp + 0.01 * exp(4.14 - 40.74 / hr) * (sysbp - diabp)
  }
}
