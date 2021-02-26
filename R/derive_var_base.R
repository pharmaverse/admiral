#' Derive BASE
#'
#' Derive the `BASE` variable in a BDS dataset
#'
#' @param bds_dataset `data.frame`
#' @param by `character` vector of columns uniqely identifying a set
#'        of records for which to calculate `BASE`
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
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
#'   "TEST01", "PAT01",  "PARAM01",  9.7,  "N",
#'   "TEST01", "PAT01",  "PARAM01", 15.01, "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
#'   "TEST01", "PAT01",  "PARAM02", NA,    "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "N",
#' )
#' derive_var_base(bds_dataset)
#'
derive_var_base <- function(bds_dataset, by = c("USUBJID", "PARAMCD", "BASETYPE")) {
  derive_baseline(bds_dataset, by, AVAL, BASE)
}

#' Derive BASEC
#'
#' Derive the `BASEC` variable in a BDS dataset
#'
#' @inheritParams derive_var_base
#'
#' @return
#' A new `data.frame` containing all records and variables of the input
#' dataset plus the `BASEC` variable
#'
#' @export
#'
#' @author Thomas Neitmann
#'
#' @examples
#' bds_dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", "LOW,     "Y",
#'   "TEST01", "PAT01",  "PARAM01", "LOW,     "N",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "N",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "N",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "N",
#' )
#' derive_var_basec(bds_dataset)
#'
derive_var_basec <- function(bds_dataset, by = c("USUBJID", "PARAMCD", "BASETYPE")) {
  derive_baseline(bds_dataset, by, AVALC, BASEC)
}

derive_baseline <- function(bds_dataset, by, source, target) {
  assert_has_variables(
    bds_dataset,
    c(by, deparse(substitute(source)), "ABLFL")
  )

  base <- bds_dataset %>%
    filter(ABLFL == "Y") %>%
    select(!!!syms(by), !!enquo(target) := !!enquo(source))

  left_join(bds_dataset, base, by = by)
}
