#' Adds a Parameter for Framingham Heart Study Cardiovascular Disease
#' 10-Year Risk Score
#'
#' Adds a record for framingham score (FCVD101) for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
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
#' @permitted [char_scalar]
#'
#' @param chol_code Total serum cholesterol code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total cholesterol assessments. This must be measured in mg/dL.
#'
#' @permitted [char_scalar]
#'
#' @param cholhdl_code HDL serum cholesterol code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HDL cholesterol assessments. This must be measured in mg/dL.
#'
#' @permitted [char_scalar]
#'
#' @param age Subject age
#'
#'   A variable containing the subject's age.
#'
#' @permitted A numeric variable name that refers to a subject age
#'                     column of the input dataset
#'
#' @param sex Subject sex
#'
#'   A variable containing the subject's sex.
#'
#' @permitted A character variable name that refers to a subject sex
#'                     column of the input dataset
#'
#' @param smokefl Smoking status flag
#'
#'   A flag indicating smoking status.
#'
#' @permitted A character variable name that refers to a smoking status
#'                     column of the input dataset.
#'
#' @param diabetfl Diabetic flag
#'
#'   A flag indicating diabetic status.
#'
#' @permitted A character variable name that refers to a diabetic
#'                     status column of the input dataset
#'
#' @param trthypfl Treated with hypertension medication flag
#'
#'   A flag indicating if a subject was treated with hypertension medication.
#'
#' @permitted A character variable name that refers to a column that
#'                     indicates whether a subject is treated for high blood
#'                     pressure
#'
#' @inheritParams derive_param_qtc
#'
#' @details
#' The values of `age`, `sex`, `smokefl`, `diabetfl` and `trthypfl` will be
#' added to the `by_vars` list.
#' The predicted probability of having cardiovascular disease (CVD)
#' within 10-years according to Framingham formula. See AHA Journal article
#' General Cardiovascular Risk Profile for Use in Primary Care for reference.
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
#' \deqn{Risk = 100 * (1 - RiskPeriodFactor^{RiskFactors})}
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords der_prm_bds_findings
#' @family der_prm_bds_findings
#'
#' @export
#'
#' @seealso [compute_framingham()]
#'
#' @examples
#' library(tibble)
#'
#' adcvrisk <- tribble(
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
#'   "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 65.55,
#'   "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
#'   "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 44.62,
#'   "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y"
#' )
#'
#'
#' adcvrisk %>%
#'   derive_param_framingham(
#'     by_vars = exprs(USUBJID, VISIT),
#'     set_values_to = exprs(
#'       PARAMCD = "FCVD101",
#'       PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
#'     ),
#'     get_unit_expr = AVALU
#'   )
#'
#' derive_param_framingham(
#'   adcvrisk,
#'   by_vars = exprs(USUBJID, VISIT),
#'   set_values_to = exprs(
#'     PARAMCD = "FCVD101",
#'     PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_framingham <- function(dataset,
                                    by_vars,
                                    set_values_to = exprs(PARAMCD = "FCVD101"),
                                    sysbp_code = "SYSBP",
                                    chol_code = "CHOL",
                                    cholhdl_code = "CHOLHDL",
                                    age = AGE,
                                    sex = SEX,
                                    smokefl = SMOKEFL,
                                    diabetfl = DIABETFL,
                                    trthypfl = TRTHYPFL,
                                    get_unit_expr,
                                    filter = NULL) {
  assert_vars(by_vars)
  age <- assert_symbol(enexpr(age))
  sex <- assert_symbol(enexpr(sex))
  smokefl <- assert_symbol(enexpr(smokefl))
  diabetfl <- assert_symbol(enexpr(diabetfl))
  trthypfl <- assert_symbol(enexpr(trthypfl))

  assert_data_frame(
    dataset,
    required_vars = expr_c(
      by_vars,
      exprs(PARAMCD, AVAL),
      age,
      sex,
      smokefl,
      diabetfl,
      trthypfl
    )
  )

  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(sysbp_code)
  assert_character_scalar(chol_code)
  assert_character_scalar(cholhdl_code)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
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
      age = !!age,
      sex = !!sex,
      smokefl = !!smokefl,
      diabetfl = !!diabetfl,
      trthypfl = !!trthypfl
    )
  )


  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(sysbp_code, chol_code, cholhdl_code),
    by_vars = expr_c(
      by_vars,
      age,
      sex,
      smokefl,
      diabetfl,
      trthypfl
    ),
    set_values_to = exprs(
      AVAL = !!analysis_value,
      !!!set_values_to
    )
  )
}
