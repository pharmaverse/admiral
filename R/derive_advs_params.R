#' Adds a parameter for BSA (Body Surface Area) using the specified method
#'
#' Adds a record for BSA (Body Surface Area) using the specified method for each by group
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
#' @param method Derivation method to use
#'
#'   The derivation method, e.g. Mosteller will use sqrt(height(cm) * weight(kg)) / 3600
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
#'   method = "Mosteller")
derive_param_bsa <- function(dataset,
                             by_vars,
                             set_values_to = vars(PARAMCD = "BSA", PARAM = "Body Surface Area", AVALU = "m^2"),
                             method = "Mosteller",
                             height_code = "HEIGHT",
                             weight_code = "WEIGHT",
                             unit_var = NULL,
                             filter = NULL) {

  assert_data_frame(dataset,
                    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL, AVALU), unit_var))
  assert_vars(by_vars)
  assert_character_scalar(method, values = c("Mosteller", "DuBois-DuBois", "Haycock",
                                             "Gehan-George", "Boyd", "Fujimoto",
                                             "Takahira"))
  assert_character_scalar(height_code)
  assert_character_scalar(weight_code)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

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

  if (method == "Mosteller") {
    bsa_formula <- expr(sqrt(!!sym(paste0("AVAL.", height_code)) *
                              !!sym(paste0("AVAL.", weight_code)) / 3600))
  }
  else if (method == "DuBois-DuBois") {
    # Note: the DuBois & DuBois formula expects the value of height in meters; we need to convert from cm.
    bsa_formula <- expr(0.20247 * (!!sym(paste0("AVAL.", height_code)) /100) ^ 0.725 *
                                   (!!sym(paste0("AVAL.", weight_code))) ^ 0.425)
  }
  else if (method == "Haycock") {
    bsa_formula <- expr(0.024265 * (!!sym(paste0("AVAL.", height_code))) ^ 0.3964 *
                                    (!!sym(paste0("AVAL.", weight_code))) ^ 0.5378)
  }
  else if (method == "Gehan-George") {
    bsa_formula <- expr(0.0235 * (!!sym(paste0("AVAL.", height_code))) ^ 0.42246 *
                                  (!!sym(paste0("AVAL.", weight_code))) ^ 0.51456)
  }
  else if (method == "Boyd") {
    # Note: the Boyd formula expects the value of weight in grams; we need to convert from kg.
    bsa_formula <- expr(0.0003207 * (!!sym(paste0("AVAL.", height_code))) ^ 0.3 *
                                     (1000 * !!sym(paste0("AVAL.", weight_code))) ^ (0.7285 - (0.0188 * log10(1000 * !!sym(paste0("AVAL.", weight_code))))))
  }
  else if (method == "Fujimoto") {
    bsa_formula <- expr(0.008883 * (!!sym(paste0("AVAL.", height_code))) ^ 0.663 *
                          (!!sym(paste0("AVAL.", weight_code))) ^ 0.444)
  }
  else if (method == "Takahira") {
    bsa_formula <- expr(0.007241 * (!!sym(paste0("AVAL.", height_code))) ^ 0.725 *
                          (!!sym(paste0("AVAL.", weight_code))) ^ 0.425)
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(height_code, weight_code),
    by_vars = by_vars,
    analysis_value = !!bsa_formula,
    set_values_to = vars(!!!set_values_to)
  )
}
