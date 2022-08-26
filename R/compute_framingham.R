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
#' compute_framingham(
#'   sysbp = 133, chol = 216.16, cholhdl = 54.91, age = 53,
#'   sex = "M", smokefl = "N", diabetfl = "N", trthypfl = "N"
#' )
#'
#' compute_framingham(
#'   sysbp = 161, chol = 186.39, cholhdl = 64.19, age = 52,
#'   sex = "F", smokefl = "Y", diabetfl = "N", trthypfl = "Y"
#' )
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
      1 - (0.95012^exp((2.32888 * log(age))
                       + (1.20904 * log(chol))
                       - (0.70833 * log(cholhdl))
                       + (2.76157 * log(if_else(trthypfl == "N", sysbp, 1)))
                       + (2.82263 * log(if_else(trthypfl == "Y", sysbp, 1)))
                       + (0.52873 * (if_else(smokefl == "Y", 1, 0)))
                       + (0.69154 * (if_else(diabetfl == "Y", 1, 0)))
                       - 26.1931)),
    sex == "M" ~
      1 - (0.88936^exp((3.06117 * log(age))
                       + (1.12370 * log(chol))
                       - (0.93263 * log(cholhdl))
                       + (1.93303 * log(if_else(trthypfl == "N", sysbp, 1)))
                       + (1.99881 * log(if_else(trthypfl == "Y", sysbp, 1)))
                       + (0.65451 * (if_else(smokefl == "Y", 1, 0)))
                       + (0.57367 * (if_else(diabetfl == "Y", 1, 0)))
                       - 23.9802))
  )


  aval * 100
}
