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
#' @param bds_dataset `data.frame`. Required variables are `AVAL` and
#' `BASE`.
#'
#' @details
#' Percent change from baseline is calculated by dividing change from
#' baseline by the absolute value of the baseline value and
#' multiplying the result by `100`.
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with an additional column named `PCHG`
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
#' derive_var_pchg(advs)
#'
derive_var_pchg <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("AVAL", "BASE"))

  bds_dataset %>%
    mutate(PCHG = if_else(BASE == 0, NA_real_, (AVAL - BASE) / abs(BASE) * 100))
}
