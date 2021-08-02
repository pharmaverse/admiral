#' Adds a parameter for BSA (Body Surface Area) using the Mosteller formula
#'
#' Adds a record for BSA (Body Surface Area) using the Mosteller formula for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived:
#' \deqn{\sqrt{HEIGHT * WEIGHT/3600}}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and the `unit_var` parameter,
#'   `PARAMCD`, and `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `HEIGHT` and `WEIGHT`.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param height_code HEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HEIGHT assessments. It is expected that HEIGHT is measured in cm.
#'
#'   Permitted Values: character value
#'
#' @param weight_code WEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the WEIGHT assessments. It is expected that WEIGHT is measured in kg.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param unit_var Variable providing the unit of the parameter
#'
#'   The variable is used to check the units of the input parameters and it is
#'   set to `"m^2"` for the new parameter.
#'
#'   Permitted Values: A variable of the input dataset
#'
#' @inheritParams derive_derived_param
#'
#' @author Eric Simms
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation advs
#'
#' @export
#'
#' @examples
#' advs <- tibble::tribble(
#' ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
#' "01-701-1015", "HEIGHT",     "Height (cm)", 170, "cm", "BASELINE",
#' "01-701-1015", "WEIGHT",     "Weight (kg)",  75, "kg", "BASELINE",
#' "01-701-1015", "WEIGHT",     "Weight (kg)",  78, "kg", "MONTH 1",
#' "01-701-1015", "WEIGHT",     "Weight (kg)",  80, "kg", "MONTH 2",
#' "01-701-1028", "HEIGHT",     "Height (cm)", 185, "cm", "BASELINE",
#' "01-701-1028", "WEIGHT",     "Weight (kg)",  90, "kg", "BASELINE",
#' "01-701-1028", "WEIGHT",     "Weight (kg)",  88, "kg", "MONTH 1",
#' "01-701-1028", "WEIGHT",     "Weight (kg)",  85, "kg", "MONTH 2",
#' )
#' derive_param_bsa(
#'   advs,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(PARAM = "Body Surface Area"))
derive_param_bsa <- function(dataset,
                             filter = NULL,
                             new_param = "BSA",
                             height_code = "HEIGHT",
                             weight_code = "WEIGHT",
                             by_vars,
                             unit_var = NULL,
                             set_values_to = vars(PARAM = "Body Surface Area", AVALU = "m^2"),
                             drop_values_from = vars(ends_with("U"))) {
  assert_character_scalar(new_param)
  assert_character_scalar(height_code)
  assert_character_scalar(weight_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset,
                    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL, AVALU), unit_var))
  assert_param_does_not_exist(dataset, new_param)
  assert_varval_list(set_values_to, optional = TRUE)

  if (!quo_is_null(unit_var)) {
    assert_unit(dataset,
                param = height_code,
                unit = "cm",
                unit_var = !!unit_var)
    assert_unit(dataset,
                param = weight_code,
                unit = "kg",
                unit_var = !!unit_var)
    set_unit_var <- vars(!!unit_var := "m^2")
  }
  else {
    set_unit_var <- NULL
  }
  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(height_code, weight_code),
    by_vars = by_vars,
    analysis_value = sqrt(!!sym(paste0("AVAL.", height_code)) * !!sym(paste0("AVAL.", weight_code)) /
                                                              3600),
    set_values_to = vars(PARAMCD = !!new_param,
                         !!!set_unit_var,
                         !!!set_values_to),
    drop_values_from = drop_values_from
  )
}
