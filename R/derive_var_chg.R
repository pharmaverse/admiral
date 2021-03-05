#' Derive Change from Baseline
#'
#' Derive change from baseline (`CHG`) in a BDS dataset
#'
#' @param bds_dataset `data.frame`. Required variables are `AVAL` and
#' `BASE`.
#'
#' @details
#' Change from baseline is calculated by subtracting the baseline value
#' from the analysis value.
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with an additional column named `CHG`
#' @export
#'
#' @examples
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
#'   "P01",    "WEIGHT", 80,    "Y",    80,
#'   "P01",    "WEIGHT", 80.8,  "",     80,
#'   "P01",    "WEIGHT", 81.4,  "",     80,
#'   "P02",    "WEIGHT", 75.3,  "Y",    75.3,
#'   "P02",    "WEIGHT", 76,    "",     75.3
#' )
#' derive_var_chg(advs)
#'
derive_var_chg <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("AVAL", "BASE"))

  bds_dataset %>%
    mutate(CHG = AVAL - BASE)
}

#' Derive Percent Change from Baseline
#'
#' Derive percent change from baseline (`PCHG`) in a BDS dataset
#'
#' @param bds_dataset `data.frame`. Required variables are `BASE` and
#' `CHG`.
#'
#' @details
#' Percent change from baseline is calculated by dividing change from
#' baseline by the baseline value.
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with an additional column named `PCHG`
#' @export
#'
#' @examples
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE, ~CHG,
#'   "P01",    "WEIGHT", 80,    "Y",    80,    0,
#'   "P01",    "WEIGHT", 80.8,  "",     80,    0.8,
#'   "P01",    "WEIGHT", 81.4,  "",     80,    1.4,
#'   "P02",    "WEIGHT", 75.3,  "Y",    75.3,  0,
#'   "P02",    "WEIGHT", 76,    "",     75.3,  0.7
#' )
#' derive_var_pchg(advs)
#'
derive_var_pchg <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("BASE", "CHG"))

  bds_dataset %>%
    mutate(PCHG = CHG / BASE)
}
