#' Derive BASE
#'
#' Derive the `BASE` variable in a BDS dataset
#'
#' @param bds_dataset `data.frame`
#'
#' @return
#' A new `data.frame` containing all records and variables of the input
#' dataset plus the `BASE` variable
#'
#' @export
#'
#' @author Thomas Neitmann
#'
#' @examples
#' bds_dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL,
#'   "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
#'   "TEST01", "PAT01",  "PARAM01",  9.7,  "N",
#'   "TEST01", "PAT01",  "PARAM01", 15.01, "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
#'   "TEST01", "PAT01",  "PARAM02", NA,    "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "N",
#' )
#' derive_var_base(bds_dataset)
#'
derive_var_base <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("STUDYID", "USUBJID", "PARAMCD", "AVAL", "ABL01FL"))

  base <- bds_dataset %>%
    filter(ABL01FL == "Y") %>%
    select(STUDYID, USUBJID, PARAMCD, BASE = AVAL)

  left_join(bds_dataset, base, by = c("STUDYID", "USUBJID", "PARAMCD"))
}
