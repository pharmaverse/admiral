#' Adds a Parameter for BSA (Body Surface Area) Using the Specified Method
#'
#' @description Adds a record for BSA (Body Surface Area) using the specified derivation
#' method for each by group (e.g., subject and visit) where the source parameters are
#' available.
#'
#' **Note:** This is a wrapper function for the more generic `derive_param_computed()`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `HEIGHT` and `WEIGHT`.
#'
#' @param method Derivation method to use. Note that `HEIGHT` is expected
#'    in cm and `WEIGHT` is expected in kg:
#'
#'   Mosteller: `sqrt(height * weight / 3600)`
#'
#'   DuBois-DuBois: `0.20247 * (height/100) ^ 0.725 * weight ^ 0.425`
#'
#'   Haycock: `0.024265 * height ^ 0.3964 * weight ^ 0.5378`
#'
#'   Gehan-George: `0.0235 * height ^ 0.42246 * weight ^ 0.51456`
#'
#'   Boyd: `0.0003207 * (height ^ 0.3) * (1000 * weight) ^
#'                  (0.7285 - (0.0188 * log10(1000 * weight)))`
#'
#'   Fujimoto: `0.008883 * height ^ 0.663 * weight ^ 0.444`
#'
#'   Takahira: `0.007241 * height ^ 0.725 * weight ^ 0.425`
#'
#' @permitted character value
#'
#' @param height_code HEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HEIGHT assessments. It is expected that HEIGHT is measured in cm.
#'
#' @permitted character value
#'
#' @param weight_code WEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the WEIGHT assessments. It is expected that WEIGHT is measured in kg.
#'
#' @permitted character value
#'
#' @param constant_by_vars By variables for when HEIGHT is constant
#'
#'   When HEIGHT is constant, the HEIGHT parameters (measured only once) are merged
#'   to the other parameters using the specified variables.
#'
#'   If height is constant (e.g. only measured once at screening or baseline) then
#'   use `constant_by_vars` to select the subject-level variable to merge on (e.g. `USUBJID`).
#'   This will produce BSA at all visits where weight is measured.  Otherwise
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
#' @seealso [compute_bsa()]
#'
#' @examples
#' library(tibble)
#'
#' # Example 1: Derive BSA where height is measured only once using constant_by_vars
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 170, "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 75, "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 78, "MONTH 1",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 80, "MONTH 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 185, "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 90, "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 88, "MONTH 1",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 85, "MONTH 2",
#' )
#'
#' derive_param_bsa(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Mosteller",
#'   set_values_to = exprs(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM),
#'   constant_by_vars = exprs(USUBJID)
#' )
#'
#' derive_param_bsa(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Fujimoto",
#'   set_values_to = exprs(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM),
#'   constant_by_vars = exprs(USUBJID)
#' )
#'
#' # Example 2: Derive BSA where height is measured only once and keep only one record
#' # where both height and weight are measured.
#'
#' derive_param_bsa(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Mosteller",
#'   set_values_to = exprs(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
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
#' derive_param_bsa(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Mosteller",
#'   set_values_to = exprs(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_bsa <- function(dataset,
                             by_vars,
                             method,
                             set_values_to = exprs(PARAMCD = "BSA"),
                             height_code = "HEIGHT",
                             weight_code = "WEIGHT",
                             get_unit_expr,
                             filter = NULL,
                             constant_by_vars = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, PARAMCD, AVAL))
  assert_character_scalar(
    method,
    values = c(
      "Mosteller", "DuBois-DuBois", "Haycock", "Gehan-George",
      "Boyd", "Fujimoto", "Takahira"
    )
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(height_code)
  assert_character_scalar(weight_code)
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_vars(constant_by_vars, optional = TRUE)

  assert_unit(
    dataset,
    param = height_code,
    required_unit = "cm",
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    param = weight_code,
    required_unit = "kg",
    get_unit_expr = !!get_unit_expr
  )

  bsa_formula <- expr(
    {{ compute_bsa }}(
      height = !!sym(paste0("AVAL.", height_code)),
      weight = !!sym(paste0("AVAL.", weight_code)),
      method = !!method
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
        AVAL = !!bsa_formula,
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

#' Compute Body Surface Area (BSA)
#'
#' Computes BSA from height and weight making use of the specified derivation method
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
#' @param method Derivation method to use:
#'
#'   Mosteller: sqrt(height * weight / 3600)
#'
#'   DuBois-DuBois: 0.007184 * height ^ 0.725 * weight ^ 0.425
#'
#'   Haycock: 0.024265 * height ^ 0.3964 * weight ^ 0.5378
#'
#'   Gehan-George: 0.0235 * height ^ 0.42246 * weight ^ 0.51456
#'
#'   Boyd: 0.0003207 * (height ^ 0.3) * (1000 * weight) ^ (0.7285 - (0.0188 * log10(1000 * weight)))
#'
#'   Fujimoto: 0.008883 * height ^ 0.663 * weight ^ 0.444
#'
#'   Takahira: 0.007241 * height ^ 0.725 * weight ^ 0.425
#'
#' @permitted character value
#'
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The BSA (Body Surface Area) in m^2.
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @seealso [derive_param_bsa()]
#'
#' @examples
#' # Derive BSA by the Mosteller method
#' compute_bsa(
#'   height = 170,
#'   weight = 75,
#'   method = "Mosteller"
#' )
#'
#' # Derive BSA by the DuBois & DuBois method
#' compute_bsa(
#'   height = c(170, 185),
#'   weight = c(75, 90),
#'   method = "DuBois-DuBois"
#' )
compute_bsa <- function(height = height,
                        weight = weight,
                        method) {
  assert_numeric_vector(height)
  assert_numeric_vector(weight)
  assert_character_scalar(
    method,
    values = c(
      "Mosteller", "DuBois-DuBois", "Haycock", "Gehan-George",
      "Boyd", "Fujimoto", "Takahira"
    )
  )

  if (method == "Mosteller") {
    bsa <- sqrt(height * weight / 3600)
  } else if (method == "DuBois-DuBois") {
    bsa <- 0.007184 * height^0.725 * weight^0.425
  } else if (method == "Haycock") {
    bsa <- 0.024265 * height^0.3964 * weight^0.5378
  } else if (method == "Gehan-George") {
    bsa <- 0.0235 * height^0.42246 * weight^0.51456
  } else if (method == "Boyd") {
    # The Boyd formula expects the value of weight in grams
    # we need to convert from kg
    bsa <- 0.0003207 * (height^0.3) *
      (1000 * weight)^(0.7285 - (0.0188 * log10(1000 * weight))) # nolint
  } else if (method == "Fujimoto") {
    bsa <- 0.008883 * height^0.663 * weight^0.444
  } else if (method == "Takahira") {
    bsa <- 0.007241 * height^0.725 * weight^0.425
  }

  bsa
}
