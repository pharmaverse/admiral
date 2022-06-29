#' Derive Percent Change from Baseline
#'
#' Derive percent change from baseline (`PCHG`) in a BDS dataset
#'
#' @param dataset The input dataset. Required variables are `AVAL` and
#' `BASE`.
#'
#' @details
#' Percent change from baseline is calculated by dividing change from
#' baseline by the absolute value of the baseline value and
#' multiplying the result by `100`.
#'
#' @return The input dataset with an additional column named `PCHG`
#'
#' @author Thomas Neitmann
#'
#' @keywords bds derivation
#'
#'  @seealso [derive_var_chg()]
#'
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
#' derive_var_pchg(advs)
derive_var_pchg <- function(dataset) {
  assert_data_frame(dataset, required_vars = vars(AVAL, BASE))

  dataset %>%
    mutate(PCHG = if_else(BASE == 0, NA_real_, (AVAL - BASE) / abs(BASE) * 100))
}
