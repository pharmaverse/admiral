#' Derive Change from Baseline
#'
#' Derive change from baseline (`CHG`) in a BDS dataset
#'
#' @param dataset The input dataset. Required variables are `AVAL` and
#' `BASE`.
#'
#' @details
#' Change from baseline is calculated by subtracting the baseline value
#' from the analysis value.
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with an additional column named `CHG`
#'
#' @family der_bds_findings
#'
#' @keywords der_bds_findings
#' @export
#'
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
#'   "P01",    "WEIGHT", 80,    "Y",    80,
#'   "P01",    "WEIGHT", 80.8,  "",     80,
#'   "P01",    "WEIGHT", 81.4,  "",     80,
#'   "P02",    "WEIGHT", 75.3,  "Y",    75.3,
#'   "P02",    "WEIGHT", 76,    "",     75.3
#' )
#' derive_var_chg(advs)
derive_var_chg <- function(dataset) {
  assert_data_frame(dataset, required_vars = vars(AVAL, BASE))

  dataset %>%
    mutate(CHG = AVAL - BASE)
}
