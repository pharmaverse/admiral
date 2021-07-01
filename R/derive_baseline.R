#' Derive BASE
#'
#' Derive the `BASE` variable in a BDS dataset
#'
#' @param dataset The input dataset
#' @param by_vars Grouping variables uniquely identifying a set
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
#' @keywords bds derivation
#'
#' @examples
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
#'   "TEST01", "PAT01",  "PARAM01",  9.7,  "N",
#'   "TEST01", "PAT01",  "PARAM01", 15.01, "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
#'   "TEST01", "PAT01",  "PARAM02", NA,    "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "N"
#' )
#' derive_var_base(dataset, by_vars = vars(USUBJID, PARAMCD))
derive_var_base <- function(dataset, by_vars) {
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
#' @keywords bds derivation
#'
#' @examples
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "N",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "N",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "N",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "N"
#' )
#' derive_var_basec(dataset, by_vars = vars(USUBJID, PARAMCD))
derive_var_basec <- function(dataset, by_vars) {
  derive_baseline(dataset, by_vars = by_vars, source_var = AVALC, new_var = BASEC)
}

#' Derive Baseline
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
#' @export
#'
#' @export
#'
#' @author Thomas Neitmann
#'
#' @keywords bds derivation
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
#'   by_vars = vars(USUBJID, PARAMCD, BASETYPE),
#'   source_var = AVALC,
#'   new_var = BASEC
#' )
derive_baseline <- function(dataset, by_vars, source_var, new_var) {

  assert_vars(by_vars)
  assert_data_frame(dataset, by_vars)
  source_var <- assert_symbol(enquo(source_var))
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))
  assert_has_variables(dataset, c(quo_text(source_var), "ABLFL"))

  base <- dataset %>%
    filter(ABLFL == "Y") %>%
    select(!!!by_vars, !!enquo(new_var) := !!enquo(source_var))

  signal_duplicate_records(base, by_vars, "Dataset contains multiple baseline records.")

  left_join(dataset, base, by = vars2chr(by_vars))
}
