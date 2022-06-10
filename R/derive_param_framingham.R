#' Adds a Parameter for Framingham Heart Study Cardiovascular Disease
#' 10-Year Risk Score
#'
#' Adds a record for framingham score (FCVD101) for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and
#'   `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `sysbp_code`, `chol_code`
#'   and `hdl_code`.
#'
#' @param sysbp_code Systolic blood pressure parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the systolic blood pressure assessments.
#'
#'   \emph{Permitted Values:} character value
#'
#' @param chol_code Total serum cholesterol code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total cholesterol assessments. This must be measured in mg/dL.
#'
#'   \emph{Permitted Values:} character value
#'
#' @param cholhdl_code HDL serum cholesterol code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HDL cholesterol assessments. This must be measured in mg/dL.
#'
#'   \emph{Permitted Values:} character value
#'
#' @param age Subject age
#'
#'   A variable containing the subject's age.
#'
#'   \emph{Permitted Values:} A numeric variable name that refers to a subject age
#'                     column of the input dataset
#'
#' @param sex Subject sex
#'
#'   A variable containing the subject's sex.
#'
#'   \emph{Permitted Values:} A character variable name that refers to a subject sex
#'                     column of the input dataset
#'
#' @param smokefl Smoking status flag
#'
#'   A flag indicating smoking status.
#'
#'   \emph{Permitted Values:} A character variable name that refers to a smoking status
#'                     column of the input dataset.
#'
#' @param diabetfl Diabetic flag
#'
#'   A flag indicating diabetic status.
#'
#'   \emph{Permitted Values:} A character variable name that refers to a diabetic
#'                     status column of the input dataset
#'
#' @param trthypfl Treated with hypertension medication flag
#'
#'   A flag indicating if a subject was treated with hypertension medication.
#'
#'   \emph{Permitted Values:} A character variable name that refers to a column that
#'                     indicates whether a subject is treated for high blood
#'                     pressure
#'
#' @inheritParams derive_derived_param
#'
#' @inheritParams derive_param_qtc
#'
#' @details
#' The values of `age`, `sex`, `smokefl`, `diabetfl` and `trthypfl` will be
#' added to the `by_vars` list.
#'
#' The predicted probability of having cardiovascular disease (CVD)
#' within 10-years according to Framingham formula
#' \href{https://www.ahajournals.org/doi/pdf/10.1161/CIRCULATIONAHA.107.699579}{D'Agostino, 2008} is:
#'
#' \strong{For Women:}
#'
#' \tabular{rr}{
#' \strong{Factor} \tab \strong{Amount} \cr
#' Age \tab 2.32888 \cr
#' Total Chol \tab 1.20904 \cr
#' HDL Chol \tab -0.70833 \cr
#' Sys BP \tab 2.76157 \cr
#' Sys BP + Hypertension Meds \tab 2.82263 \cr
#' Smoker \tab 0.52873 \cr
#' Non-Smoker \tab 0 \cr
#' Diabetic \tab 0.69154 \cr
#' Not Diabetic \tab 0 \cr
#' Average Risk \tab 26.1931 \cr
#' Risk Period \tab 0.95012 \cr
#' }
#'
#'\strong{For Men:}
#'
#' \tabular{rr}{
#' \strong{Factor} \tab \strong{Amount} \cr
#' Age \tab 3.06117 \cr
#' Total Chol \tab 1.12370 \cr
#' HDL Chol \tab -0.93263 \cr
#' Sys BP \tab 1.93303 \cr
#' Sys BP + Hypertension Meds \tab 2.99881 \cr
#' Smoker \tab .65451  \cr
#' Non-Smoker \tab 0 \cr
#' Diabetic \tab 0.57367  \cr
#' Not Diabetic \tab 0 \cr
#' Average Risk \tab 23.9802 \cr
#' Risk Period \tab 0.88936 \cr
#' }
#'
#' \strong{The equation for calculating risk:}
#'
#' \deqn{RiskFactors = (log(Age) * AgeFactor)
#' + (log(TotalChol) * TotalCholFactor)
#' + (log(CholHDL) * CholHDLFactor) \\
#' + (log(SysBP) * SysBPFactor) + Smoker
#' + Diabetes Present - AvgRisk}
#'
#' \deqn{Risk = 100 * (1 - RiskPeriodFactor ^ exp(RiskFactors))}
#'
#' @author
#' Alice Ehmann
#' Jack McGavigan
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation bds
#'
#' @export
#'
#' @seealso [compute_framingham()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adcvrisk <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU,
#'   ~VISIT, ~AGE, ~SEX, ~SMOKEFL, ~DIABETFL, ~TRTHYPFL,
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,
#'   "mmHg", "BASELINE", 44, "F", "N", "N", "N",
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 115,
#'   "mmHg", "WEEK 2", 44, "F", "N", "N", "Y",
#'   "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 216.16,
#'   "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
#'   "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 210.78,
#'   "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
#'   "01-701-1015", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 54.91,
#'   "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
#'   "01-701-1015", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 26.72,
#'   "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 119,
#'   "mmHg", "BASELINE", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 101,
#'   "mmHg", "WEEK 2", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 292.01,
#'   "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 246.73,
#'   "mg/dL", "WEEK 2", 55,  "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 65.55,
#'   "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 44.62,
#'   "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y"
#' )
#'
#'
#' adcvrisk %>%
#'   derive_param_framingham(
#'     by_vars = vars(USUBJID, VISIT),
#'     set_values_to = vars(
#'       PARAMCD = "FCVD101",
#'       PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
#'     ),
#'     get_unit_expr = AVALU
#'   )
#'
#' derive_param_framingham(
#'   adcvrisk,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "FCVD101",
#'     PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#'

derive_param_framingham <- function(dataset,
                                    by_vars,
                                    set_values_to = vars(PARAMCD = "FCVD101"),
                                    sysbp_code = "SYSBP",
                                    chol_code = "CHOL",
                                    cholhdl_code = "CHOLHDL",
                                    age = AGEdev,
                                    sex = SEX,
                                    smokefl = SMOKEFL,
                                    diabetfl = DIABETFL,
                                    trthypfl = TRTHYPFL,
                                    get_unit_expr,
                                    filter = NULL) {

  assert_symbol(enquo(AGE))
  assert_symbol(enquo(SEX))
  assert_symbol(enquo(SMOKEFL))
  assert_symbol(enquo(DIABETFL))
  assert_symbol(enquo(TRTHYPFL))

  assert_vars(by_vars)
  assert_data_frame(dataset,
                    required_vars = quo_c(vars(!!!by_vars), PARAMCD, AVAL))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(sysbp_code)
  assert_character_scalar(chol_code)
  assert_character_scalar(cholhdl_code)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  assert_unit(
    dataset,
    sysbp_code,
    required_unit = "mmHg",
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    chol_code,
    required_unit = "mg/dL",
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    cholhdl_code,
    required_unit = "mg/dL",
    get_unit_expr = !!get_unit_expr
  )

  analysis_value <- expr(
    compute_framingham(
      sysbp = !!sym(paste0("AVAL.", sysbp_code)),
      chol = !!sym(paste0("AVAL.", chol_code)),
      cholhdl = !!sym(paste0("AVAL.", cholhdl_code)),
      age = !!enquo(age),
      sex = !!enquo(sex),
      smokefl = !!enquo(smokefl),
      diabetfl = !!enquo(diabetfl),
      trthypfl = !!enquo(trthypfl)
    )
  )

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(sysbp_code, chol_code, cholhdl_code),
    by_vars = quo_c(
      vars(!!!by_vars),
      enquo(age),
      enquo(sex),
      enquo(smokefl),
      enquo(diabetfl),
      enquo(trthypfl)
    ),
    analysis_value = !!analysis_value,
    set_values_to = set_values_to
  )

}

#' Compute Framingham Heart Study Cardiovascular Disease 10-Year Risk Score
#'
#' Computes Framingham Heart Study Cardiovascular Disease 10-Year Risk Score
#' (FCVD101) based on systolic blood pressure, total serum cholesterol (mg/dL),
#' HDL serum cholesterol (mg/dL), sex, smoking status, diabetic status,
#' and treated for hypertension flag.
#'
#' @param sysbp Systolic blood pressure
#'
#'   A numeric vector is expected.
#'
#' @param chol Total serum cholesterol (mg/dL)
#'
#'   A numeric vector is expected.
#'
#' @param cholhdl HDL serum cholesterol (mg/dL)
#'
#'   A numeric vector is expected.
#'
#' @param age Age (years)
#'
#'   A numeric vector is expected.
#'
#' @param sex Gender
#'
#'   A character vector is expected.
#'   Expected Values: 'M' 'F'
#'
#' @param smokefl Smoking Status
#'
#'   A character vector  is expected.
#'   Expected Values: 'Y' 'N'
#'
#' @param diabetfl Diabetic Status
#'
#'   A character vector is expected.
#'   Expected Values: 'Y' 'N'
#'
#' @param trthypfl Treated for hypertension status
#'
#'   A character vector is expected.
#'   Expected Values: 'Y' 'N'
#'
#' @author Alice Ehmann
#'
#' @details
#' The predicted probability of having cardiovascular disease (CVD)
#' within 10-years according to Framingham formula
#' \href{https://www.ahajournals.org/doi/pdf/10.1161/CIRCULATIONAHA.107.699579}{D'Agostino, 2008} is:
#'
#' \strong{For Women:}
#'
#' \tabular{rr}{
#' \strong{Factor} \tab \strong{Amount} \cr
#' Age \tab 2.32888 \cr
#' Total Chol \tab 1.20904 \cr
#' HDL Chol \tab -0.70833 \cr
#' Sys BP \tab 2.76157 \cr
#' Sys BP + Hypertension Meds \tab 2.82263 \cr
#' Smoker \tab 0.52873 \cr
#' Non-Smoker \tab 0 \cr
#' Diabetic \tab 0.69154 \cr
#' Not Diabetic \tab 0 \cr
#' Average Risk \tab 26.1931 \cr
#' Risk Period \tab 0.95012 \cr
#' }
#'
#' \strong{For Men:}
#'
#' \tabular{rr}{
#' \strong{Factor} \tab \strong{Amount} \cr
#' Age \tab 3.06117 \cr
#' Total Chol \tab 1.12370 \cr
#' HDL Chol \tab -0.93263 \cr
#' Sys BP \tab 1.93303 \cr
#' Sys BP + Hypertension Meds \tab 2.99881 \cr
#' Smoker \tab .65451  \cr
#' Non-Smoker \tab 0 \cr
#' Diabetic \tab 0.57367  \cr
#' Not Diabetic \tab 0 \cr
#' Average Risk \tab 23.9802 \cr
#' Risk Period \tab 0.88936 \cr
#' }
#'
#' \strong{The equation for calculating risk:}
#'
#' \deqn{RiskFactors = (log(Age) * AgeFactor)
#' + (log(TotalChol) * TotalCholFactor)
#' + (log(CholHDL) * CholHDLFactor) \\
#' + (log(SysBP) * SysBPFactor) + Smoker
#' + Diabetes Present - AvgRisk}
#'
#' \deqn{Risk = 100 * (1 - RiskPeriodFactor ^ exp(RiskFactors))}
#'
#' @return A numeric vector of Framingham values
#'
#' @keywords computation adam
#'
#' @export
#'
#' @seealso [derive_param_framingham()]
#'
#' @examples
#' compute_framingham(sysbp = 133, chol = 216.16, cholhdl = 54.91, age=53,
#'                    sex = "M", smokefl = "N", diabetfl = "N", trthypfl = "N")
#'
#' compute_framingham(sysbp = 161, chol = 186.39, cholhdl = 64.19, age=52,
#'                    sex = "F", smokefl = "Y", diabetfl = "N", trthypfl = "Y")
#'
compute_framingham <- function(sysbp, chol, cholhdl, age, sex, smokefl,
                               diabetfl, trthypfl) {
  assert_numeric_vector(sysbp)
  assert_numeric_vector(chol)
  assert_numeric_vector(cholhdl)
  assert_numeric_vector(age)
  assert_character_vector(sex, values = c("M", "F"))
  assert_character_vector(smokefl, values = c("Y", "N"))
  assert_character_vector(diabetfl, values = c("Y", "N"))
  assert_character_vector(trthypfl, values = c("Y", "N"))

  aval <- case_when(
    sex == "F" ~
      1 - (0.95012 ^ exp((2.32888 * log(age))
                         + (1.20904 * log(chol))
                         - (0.70833 * log(cholhdl))
                         + (2.76157 * log(if_else(trthypfl == "N", sysbp, 1)))
                         + (2.82263 * log(if_else(trthypfl == "Y", sysbp, 1)))
                         + (0.52873 * (if_else(smokefl == "Y", 1, 0)))
                         + (0.69154 * (if_else(diabetfl == "Y", 1, 0)))
                         - 26.1931)),
    sex == "M" ~
      1 - (0.88936 ^ exp((3.06117 * log(age))
                         + (1.12370 * log(chol))
                         - (0.93263 * log(cholhdl))
                         + (1.93303 * log(if_else(trthypfl == "N", sysbp, 1)))
                         + (1.99881 * log(if_else(trthypfl == "Y", sysbp, 1)))
                         + (0.65451 * (if_else(smokefl == "Y", 1, 0)))
                         + (0.57367 * (if_else(diabetfl == "Y", 1, 0)))
                         - 23.9802)))


 aval * 100

}

