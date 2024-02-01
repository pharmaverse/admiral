#' Compute Estimated Glomerular Filtration Rate (eGFR) for Kidney Function
#'
#' Compute Kidney Function Tests:
#'   - Estimated Creatinine Clearance (CRCL) by Cockcroft-Gault equation
#'   - Estimated Glomerular Filtration Rate (eGFR) by CKD-EPI or MDRD equations
#'
#' @param creat Creatinine
#'
#'   A numeric vector is expected.
#'
#' @param creatu Creatinine Units
#'
#'   A character vector is expected.
#'
#'   Default: `"SI"`
#'
#'   Expected Values: `"SI"`, `"CV"`, `"umol/L"`, `"mg/dL"`
#'
#' @param age Age (years)
#'
#'   A numeric vector is expected.
#'
#' @param weight Weight (kg)
#'
#'   A numeric vector is expected if `method = "CRCL"`
#'
#' @param sex Gender
#'
#'   A character vector is expected.
#'
#'   Expected Values: `"M"`, `"F"`
#'
#' @param race Race
#'
#'   A character vector is expected if `method = "MDRD"`
#'
#'   Expected Values: `"BLACK OR AFRICAN AMERICAN"` and others
#'
#' @param method Method
#'
#'   A character vector is expected.
#'
#'   Expected Values: `"CRCL"`, `"CKD-EPI"`, `"MDRD"`
#'
#' @details
#'
#' Calculates an estimate of Glomerular Filtration Rate (eGFR)
#'
#' \strong{CRCL Creatinine Clearance (Cockcroft-Gault)}
#'
#' For Creatinine in umol/L:
#'
#' \deqn{\frac{(140 - age) \times weight(kg) \times constant}{Serum\:Creatinine(\mu mol/L)}}
#'
#' \deqn{Constant = 1.04\:for\:females, 1.23\:for\:males}
#'
#' For Creatinine in mg/dL:
#'
#' \deqn{\frac{(140 - age) \times weight(kg) \times (0.85\:if\:female)}{72 \times
#' Serum\:Creatinine(mg/dL)}}
#'
#' units = mL/min
#'
#' \strong{CKD-EPI Chronic Kidney Disease Epidemiology Collaboration formula}
#'
#' \deqn{eGFR = 142 \times min(SCr/{\kappa}, 1)^{\alpha} \times max(SCr/{\kappa}, 1)^{-1.200}
#' \times 0.9938^{Age} \times 1.012 [if\:female]}
#'
#' SCr = standardized serum creatinine in mg/dL
#' (Note SCr(mg/dL) = Creat(umol/L) / 88.42)
#'
#' \deqn{\kappa} = 0.7 (females) or 0.9 (males)
#' \deqn{\alpha} = -0.241 (female) or -0.302 (male)
#' units = mL/min/1.73 m2
#'
#' \strong{MDRD Modification of Diet in Renal Disease formula}
#'
#' \deqn{eGFR = 175 \times (SCr)^{-1.154} \times (age)^{-0.203}
#' \times 0.742 [if\:female] \times 1.212 [if\:Black]}
#'
#' SCr = standardized serum creatinine in mg/dL
#' (Note SCr(mg/dL) = Creat(umol/L) / 88.42)
#'
#' units = mL/min/1.73 m2
#'
#' @return A numeric vector of egfr values
#'
#' @keywords com_bds_findings
#' @family com_bds_findings
#'
#' @export
#'
#' @examples
#' compute_egfr(
#'   creat = 90, creatu = "umol/L", age = 53, weight = 85, sex = "M", method = "CRCL"
#' )
#'
#' compute_egfr(
#'   creat = 90, creatu = "umol/L", age = 53, sex = "M", race = "ASIAN", method = "MDRD"
#' )
#'
#' compute_egfr(
#'   creat = 70, creatu = "umol/L", age = 52, sex = "F", race = "BLACK OR AFRICAN AMERICAN",
#'   method = "MDRD"
#' )
#'
#' compute_egfr(
#'   creat = 90, creatu = "umol/L", age = 53, sex = "M", method = "CKD-EPI"
#' )
#'
#'
#' base <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~AGE, ~SEX, ~RACE, ~WTBL, ~CREATBL, ~CREATBLU,
#'   "P01", "P01-1001", 55, "M", "WHITE", 90.7, 96.3, "umol/L",
#'   "P01", "P01-1002", 52, "F", "BLACK OR AFRICAN AMERICAN", 68.5, 70, "umol/L",
#'   "P01", "P01-1003", 67, "M", "BLACK OR AFRICAN AMERICAN", 85.0, 77, "umol/L",
#'   "P01", "P01-1004", 76, "F", "ASIAN", 60.7, 65, "umol/L",
#' )
#'
#' base %>%
#'   dplyr::mutate(
#'     CRCL_CG = compute_egfr(
#'       creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
#'       method = "CRCL"
#'     ),
#'     EGFR_EPI = compute_egfr(
#'       creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
#'       method = "CKD-EPI"
#'     ),
#'     EGFR_MDRD = compute_egfr(
#'       creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
#'       race = RACE, method = "MDRD"
#'     ),
#'   )
compute_egfr <- function(creat, creatu = "SI", age, weight, sex, race = NULL, method) {
  assert_numeric_vector(creat)
  assert_character_vector(creatu, values = c(
    "SI", "CV", "mg/dL", "umol/L", NA_character_,
    optional = TRUE
  ))
  assert_numeric_vector(age)
  assert_character_vector(sex, values = c("M", "F"))
  assert_character_vector(race, optional = TRUE)
  assert_character_scalar(
    method,
    values = c(
      "CRCL", "MDRD", "CKD-EPI"
    )
  )

  scr <- case_when(
    tolower(creatu) %in% c("cv", "mg/dl") ~ creat,
    TRUE ~ creat / 88.42
  )

  if (method == "MDRD") {
    assert_character_vector(race)

    egfr <- case_when(
      race == "BLACK OR AFRICAN AMERICAN" & sex == "F" ~ 175 * (scr^-1.154) *
        (age^-0.203) * 0.742 * 1.212,
      race == "BLACK OR AFRICAN AMERICAN" ~ 175 * (scr^-1.154) * (age^-0.203) * 1.212,
      sex == "F" ~ 175 * (scr^-1.154) * (age^-0.203) * 0.742,
      sex == "M" ~ 175 * (scr^-1.154) * (age^-0.203)
    )
  } else if (method == "CRCL") {
    assert_numeric_vector(weight)

    egfr <- case_when(
      tolower(creatu) %in% c("cv", "mg/dl") & sex == "F" ~
        ((140 - age) * weight * 0.85) / (creat * 72),
      tolower(creatu) %in% c("cv", "mg/dl") & sex == "M" ~
        ((140 - age) * weight) / (creat * 72),
      sex == "F" ~ ((140 - age) * weight * 1.04) / creat,
      sex == "M" ~ ((140 - age) * weight * 1.23) / creat
    )
  } else if (method == "CKD-EPI") {
    kappa <- case_when(
      sex == "F" ~ 0.7,
      sex == "M" ~ 0.9,
      TRUE ~ NA_real_
    )

    alpha <- case_when(
      sex == "F" ~ -0.241,
      sex == "M" ~ -0.302,
      TRUE ~ NA_real_
    )

    gender_coefficent <- case_when(
      sex == "F" ~ 1.012,
      TRUE ~ 1
    )

    egfr <- 142 * pmin(scr / kappa, 1)^(alpha) *
      pmax(scr / kappa, 1)^(-1.200) *
      0.9938^age *
      gender_coefficent
  }

  return(egfr)
}
