#' Adds a Parameter for Mean Arterial Pressure
#'
#' Adds a record for mean arterial pressure (MAP) for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and
#'   `AVAL` are expected.
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
#' @author Stefan Bundfuss
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
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' advs <- tribble(
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
#'     by_vars = vars(USUBJID, VISIT),
#'     set_values_to = vars(
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
#'   by_vars = vars(USUBJID, VISIT),
#'   hr_code = "PULSE",
#'   set_values_to = vars(
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_map <- function(dataset,
                             by_vars,
                             set_values_to = vars(PARAMCD = "MAP"),
                             sysbp_code = "SYSBP",
                             diabp_code = "DIABP",
                             hr_code = NULL,
                             get_unit_expr,
                             filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, PARAMCD, AVAL))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(sysbp_code)
  assert_character_scalar(diabp_code)
  assert_character_scalar(hr_code, optional = TRUE)
  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

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

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(sysbp_code, diabp_code, hr_code),
    by_vars = by_vars,
    analysis_value = !!analysis_value,
    set_values_to = set_values_to
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
#' @author Stefan Bundfuss
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

#' Adds a Parameter for BSA (Body Surface Area) Using the Specified Method
#'
#' Adds a record for BSA (Body Surface Area) using the specified derivation method
#' for each by group (e.g., subject and visit) where the source parameters are available.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and
#'   `AVAL` are expected.
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
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#' @author Eric Simms
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
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
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
#'   by_vars = vars(USUBJID, VISIT),
#'   method = "Mosteller",
#'   set_values_to = vars(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#' derive_param_bsa(
#'   advs,
#'   by_vars = vars(USUBJID, VISIT),
#'   method = "Fujimoto",
#'   set_values_to = vars(
#'     PARAMCD = "BSA",
#'     PARAM = "Body Surface Area (m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_bsa <- function(dataset,
                             by_vars,
                             method,
                             set_values_to = vars(PARAMCD = "BSA"),
                             height_code = "HEIGHT",
                             weight_code = "WEIGHT",
                             get_unit_expr,
                             filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, PARAMCD, AVAL))
  assert_character_scalar(
    method,
    values = c(
      "Mosteller", "DuBois-DuBois", "Haycock", "Gehan-George",
      "Boyd", "Fujimoto", "Takahira"
    )
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(height_code)
  assert_character_scalar(weight_code)
  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

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
    compute_bsa(
      height = !!sym(paste0("AVAL.", height_code)),
      weight = !!sym(paste0("AVAL.", weight_code)),
      method = method
    )
  )

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(height_code, weight_code),
    by_vars = by_vars,
    analysis_value = !!bsa_formula,
    set_values_to = set_values_to
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
#'   Permitted Values: numeric vector
#'
#' @param weight WEIGHT value
#'
#'   It is expected that WEIGHT is in kg.
#'
#'   Permitted Values: numeric vector
#'
#' @param method Derivation method to use:
#'
#'   Mosteller: sqrt(height * weight / 3600)
#'
#'   DuBois-DuBois: 0.20247 * (height/100) ^ 0.725 * weight ^ 0.425
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
#'   Permitted Values: character value
#'
#' @author Eric Simms
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
    # The DuBois & DuBois formula expects the value of height in meters
    # We need to convert from cm
    bsa <- 0.20247 * (height / 100)^0.725 * weight^0.425
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

#' Adds a Parameter for BMI
#'
#' Adds a record for BMI/Body Mass Index using Weight and Height each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and
#'   `AVAL` are expected.
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
#'   Permitted Values: character value
#'
#' @param height_code HEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HEIGHT. It is expected that HEIGHT is measured in cm
#'
#'   Permitted Values: character value
#'
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#' @details
#' The analysis value of the new parameter is derived as
#' \deqn{BMI = \frac{WEIGHT}{HEIGHT^2}}
#'
#' @author Pavan Kumar
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
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 147,   "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.0,  "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.4,  "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 53.1,  "WEEK 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 163,   "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 78.5,  "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.3,  "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.7,  "WEEK 2"
#' )
#'
#' derive_param_bmi(
#'   advs,
#'   by_vars = vars(USUBJID, AVISIT),
#'   weight_code = "WEIGHT",
#'   height_code = "HEIGHT",
#'   set_values_to = vars(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_bmi <- function(dataset,
                             by_vars,
                             set_values_to = vars(PARAMCD = "BMI"),
                             weight_code = "WEIGHT",
                             height_code = "HEIGHT",
                             get_unit_expr,
                             filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, PARAMCD, AVAL))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(weight_code)
  assert_character_scalar(height_code)
  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

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

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(weight_code, height_code),
    by_vars = by_vars,
    analysis_value = compute_bmi(
      height = !!sym(paste0("AVAL.", height_code)),
      weight = !!sym(paste0("AVAL.", weight_code))
    ),
    set_values_to = set_values_to
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
#'   Permitted Values: numeric vector
#'
#' @param weight WEIGHT value
#'
#'   It is expected that WEIGHT is in kg.
#'
#'   Permitted Values: numeric vector
#'
#' @author Pavan Kumar
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
#' @examples
#' compute_bmi(height = 170, weight = 75)
compute_bmi <- function(height, weight) {
  assert_numeric_vector(height)
  assert_numeric_vector(weight)

  if_else(height == 0, NA_real_, weight / (height * height / 10000))
}
