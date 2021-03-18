#' Derive BASE
#'
#' Derive the `BASE` variable in a BDS dataset
#'
#' @param dataset `data.frame`
#' @param by_vars `character` vector of columns uniqely identifying a set
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
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASETYPE,
#'   "TEST01", "PAT01",  "PARAM01", 10.12, "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01",  9.7,  "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01", 15.01, "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", NA,    "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "N",    "LAST"
#' )
#' derive_var_base(dataset)
#'
derive_var_base <- function(dataset, by_vars = c("USUBJID", "PARAMCD", "BASETYPE")) {
  derive_baseline(dataset, by_vars = by_vars, source_var = AVAL, new_var = BASE)
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
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL, ~BASETYPE,
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "N",    "LAST"
#' )
#' derive_var_basec(dataset)
#'
derive_var_basec <- function(dataset, by_vars = c("USUBJID", "PARAMCD", "BASETYPE")) {
  derive_baseline(dataset, by_vars = by_vars, source_var = AVALC, new_var = BASEC)
}

#' Derive Basline
#'
#' Derive a baseline variable, e.g. `BASE`, in a BDS dataset
#'
#' @inheritParams derive_var_base
#' @param source_var The column from which to extract the baseline value,
#' e.g. `AVAL`
#' @param new_var The name of the newly created baseline column, e.g.
#' `BASE`
#'
#' @return
#' A new `data.frame` containing all records and variables of the input
#' dataset plus the `new_var` variable.
#'
#' @author Thomas Neitmann
#'
#' @examples
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL, ~BASETYPE,
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "N",    "LAST",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "N",    "LAST"
#' )
#' derive_baseline(
#'   dataset,
#'   by_vars = c("USUBJID", "PARAMCD", "BASETYPE"),
#'   source_var = AVALC,
#'   new_var = BASEC
#' )
#'
derive_baseline <- function(dataset, by_vars, source_var, new_var) {
  assert_has_variables(
    dataset,
    c(by_vars, deparse(substitute(source_var)), "ABLFL")
  )

  base <- dataset %>%
    filter(ABLFL == "Y") %>%
    select(!!!syms(by_vars), !!enquo(new_var) := !!enquo(source_var))

  assert_has_only_one_baseline_record(base, by_vars)

  left_join(dataset, base, by = by_vars)
}
