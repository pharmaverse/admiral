#' Compute Creatinine Clearance (CRCL)
#'
#' Compute Creatinine Clearance (CRCL) as estimate of renal function
#'
#' @param creat Creatinine (umol/L)
#'
#'   A numeric vector is expected.
#'
#' @param age Age (years)
#'
#'   A numeric vector is expected.
#'
#' @param wt Weight (kg)
#'
#'   A numeric vector is expected.
#'
#' @param sex Gender
#'
#'   A character vector is expected.
#'   Expected Values: 'M' 'F'
#'
#' @details
#' Calculates an estimate of Creatinine Clearance
#'
#' \deqn{((140 - age) * weight (kg) * constant) / creatinine (umol/L)}
#'
#' Constant = 1.04 for females, 1.23 for males
#'
#' @return A numeric vector of crcl values (mL/min)
#'
#' @keywords com_bds_findings
#' @family com_bds_findings
#'
#' @export
#'
#' @examples
#' compute_crcl(
#'   creat = 90, age = 53, wt = 85, sex = "M"
#' )
#'
#' compute_crcl(
#'   creat = 70, age = 52, wt = 68, sex = "F"
#' )
compute_crcl <- function(creat, age, wt, sex) {
  assert_numeric_vector(creat)
  assert_numeric_vector(age)
  assert_numeric_vector(wt)
  assert_character_vector(sex, values = c("M", "F"))


  crcl <- if_else(sex == "F",
    ((140 - age) * wt * 1.04) / creat,
    ((140 - age) * wt * 1.23) / creat
  )
}
