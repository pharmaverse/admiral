#' Adds a Parameter for BMI
#'
#' @description Adds a record for BMI/Body Mass Index using Weight and Height each by group
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
#'   parameter) and to the parameters specified by `weight_code` and `height_code`.
#'
#' @param weight_code WEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the WEIGHT. It is expected that WEIGHT is measured in kg
#'
#' @permitted character value
#'
#' @param height_code HEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HEIGHT. It is expected that HEIGHT is measured in cm
#'
#' @permitted character value
#'
#' @permitted logical scalar
#'
#' @param constant_by_vars By variables for when HEIGHT is constant
#'
#'   When HEIGHT is constant, the HEIGHT parameters (measured only once) are merged
#'   to the other parameters using the specified variables.
#'
#'   If height is constant (e.g. only measured once at screening or baseline) then
#'   use `constant_by_vars` to select the subject-level variable to merge on (e.g. `USUBJID`).
#'   This will produce BMI at all visits where weight is measured.  Otherwise
#'   it will only be calculated at visits with both height and weight collected.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @inheritParams derive_param_map
#'
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#' @details
#' The analysis value of the new parameter is derived as
#' \deqn{BMI = \frac{WEIGHT}{HEIGHT^2}}
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
#' @seealso [compute_bmi()]
#'
#' @examples
#'
#' # Example 1: Derive BMI where height is measured only once using constant_by_vars
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 147, "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.0, "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.4, "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 53.1, "WEEK 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 163, "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 78.5, "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.3, "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.7, "WEEK 2"
#' )
#'
#' derive_param_bmi(
#'   advs,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   weight_code = "WEIGHT",
#'   height_code = "HEIGHT",
#'   set_values_to = exprs(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM),
#'   constant_by_vars = exprs(USUBJID)
#' )
#'
#' # Example 2: Derive BMI where height is measured only once and keep only one record
#' # where both height and weight are measured.
#' derive_param_bmi(
#'   advs,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   weight_code = "WEIGHT",
#'   height_code = "HEIGHT",
#'   set_values_to = exprs(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#' # Example 3: Pediatric study where height and weight are measured multiple times
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
#'   "01-101-1001", "HEIGHT", "Height (cm)", 47.1, "BASELINE",
#'   "01-101-1001", "HEIGHT", "Height (cm)", 59.1, "WEEK 12",
#'   "01-101-1001", "HEIGHT", "Height (cm)", 64.7, "WEEK 24",
#'   "01-101-1001", "HEIGHT", "Height (cm)", 68.2, "WEEK 48",
#'   "01-101-1001", "WEIGHT", "Weight (kg)", 2.6, "BASELINE",
#'   "01-101-1001", "WEIGHT", "Weight (kg)", 5.3, "WEEK 12",
#'   "01-101-1001", "WEIGHT", "Weight (kg)", 6.7, "WEEK 24",
#'   "01-101-1001", "WEIGHT", "Weight (kg)", 7.4, "WEEK 48",
#' )
#'
#' derive_param_bmi(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   weight_code = "WEIGHT",
#'   height_code = "HEIGHT",
#'   set_values_to = exprs(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_bmi <- function(dataset,
                             by_vars,
                             set_values_to = exprs(PARAMCD = "BMI"),
                             weight_code = "WEIGHT",
                             height_code = "HEIGHT",
                             get_unit_expr,
                             filter = NULL,
                             constant_by_vars = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, PARAMCD, AVAL))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(weight_code)
  assert_character_scalar(height_code)
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_vars(constant_by_vars, optional = TRUE)


  assert_unit(
    dataset,
    param = weight_code,
    required_unit = "kg",
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    param = height_code,
    required_unit = "cm",
    get_unit_expr = !!get_unit_expr
  )

  bmi_formula <- expr(
    compute_bmi(
      height = !!sym(paste0("AVAL.", height_code)),
      weight = !!sym(paste0("AVAL.", weight_code))
    )
  )

  if (is.null(constant_by_vars)) {
    parameters <- c(weight_code, height_code)
    constant_parameters <- NULL
  } else {
    parameters <- c(weight_code)
    constant_parameters <- c(height_code)
  }

  withCallingHandlers(
    derive_param_computed(
      dataset,
      filter = !!filter,
      parameters = parameters,
      by_vars = by_vars,
      set_values_to = exprs(
        AVAL = !!bmi_formula,
        !!!set_values_to
      ),
      constant_parameters = constant_parameters,
      constant_by_vars = constant_by_vars
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

#' Compute Body Mass Index (BMI)
#'
#' Computes BMI from height and weight
#'
#' @param height HEIGHT value
#'
#'   It is expected that HEIGHT is in cm.
#'
#' @permitted numeric vector
#'
#' @param weight WEIGHT value
#'
#'   It is expected that WEIGHT is in kg.
#'
#' @permitted numeric vector
#'
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The BMI (Body Mass Index Area) in kg/m^2.
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @seealso [derive_param_bmi()]
#'
#' @examples
#' compute_bmi(height = 170, weight = 75)
compute_bmi <- function(height, weight) {
  assert_numeric_vector(height)
  assert_numeric_vector(weight)

  if_else(height == 0, NA_real_, weight / (height * height / 10000))
}
