#' Adds a Parameter for BMI
#'
#' Adds a record for BMI/Body Mass Index using Weight and Height each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{BMI = {WEIGHT / (HEIGHT)^2}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and the `unit_var` parameter,
#'   `PARAMCD`, and `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `weight_code` and `height_code`.
#'
#' @param weight_code Weight parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments. It is expected that QT is measured in msec.
#'
#'   Permitted Values: character value
#'
#' @param height_code Height parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments. It is expected that RR is measured in msec.
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
#'   set to `"kg/m^2"` for the new parameter.
#'
#'   Permitted Values: A variable of the input dataset
#'
#' @inheritParams derive_derived_param
#'
#' @author Pavan Kumar
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation advs
#'
#' @export
#'
#' @examples
#'
#' # derive BMI where height is measured only once
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 147, "cm", "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.0, "kg", "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.4, "kg", "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 53.1, "kg", "WEEK 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 163, "cm", "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 78.5, "kg", "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.3, "kg", "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.7, "kg", "WEEK 2"
#' )
#'
#' derive_advs_bmi(
#'   advs,
#'   by_vars = vars(USUBJID, VISIT),
#'   weight_code = "WEIGHT",
#'   height_code = "HEIGHT",
#'   set_values_to = vars(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)",
#'     AVALU = "kg/m^2"
#'   )
#' )
#'


derive_param_bmi <-  function(dataset,
                              by_vars,
                              set_values_to = vars(PARAMCD = "BMI"),
                              weight_code = "WEIGHT",
                              height_code = "HEIGHT",
                              unit_var = NULL,
                              filter = NULL) {
  assert_character_scalar(weight_code)
  assert_character_scalar(height_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  if (!quo_is_null(unit_var)) {
    assert_unit(
      dataset,
      param = weight_code,
      unit = "kg",
      unit_var = !!unit_var
    )
    assert_unit(
      dataset,
      param = height_code,
      unit = "m",
      unit_var = !!unit_var
    )
    set_unit_var <- vars(!!unit_var := "kg/m^2")
  } else {
    set_unit_var <- NULL
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(weight_code, height_code),
    by_vars = by_vars,
    analysis_value = !!sym(paste0("AVAL.", weight_code)) / !!sym(paste0("AVAL.", height_code)) ^2,
    set_values_to = vars(!!!set_unit_var, !!!set_values_to)
  )
}
