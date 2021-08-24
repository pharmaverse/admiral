#' Derive BASE Variable
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
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~AVISIT,    ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", 10.12, "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM01",  9.7,  "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM01", 15.01, "Day 14",   "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM02", NA,    "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM02",  8.35, "Day 14",   "N"
#' )
#' derive_var_base(dataset, by_vars = vars(USUBJID, PARAMCD))
derive_var_base <- function(dataset, by_vars) {
  derive_baseline(dataset, by_vars = by_vars, source_var = AVAL, new_var = BASE)
}

#' Derive BASEC Variable
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
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~AVISIT,    ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "Day 14",   "N",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "Day 14",   "N"
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
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~AVISIT,    ~ABLFL,
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM01", "LOW",    "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM01", "MEDIUM", "Day 14",   "N",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Baseline", "Y",
#'   "TEST01", "PAT01",  "PARAM02", "HIGH",   "Day 7",    "N",
#'   "TEST01", "PAT01",  "PARAM02", "MEDIUM", "Day 14",   "N",
#'   "TEST01", "PAT02",  "PARAM01", "MEDIUM", "Baseline", "Y",
#'   "TEST01", "PAT02",  "PARAM01", "LOW",    "Day 7",    "N",
#'   "TEST01", "PAT02",  "PARAM01", "MEDIUM", "Day 14",   "N",
#' )
#'
#' derive_baseline(
#'   dataset,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   source_var = AVALC,
#'   new_var = BASEC
#' )
derive_baseline <- function(dataset, by_vars, source_var, new_var) {
  by_vars <- assert_vars(by_vars)
  source_var <- assert_symbol(enquo(source_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, source_var, quo(ABLFL))
  )
  warn_if_vars_exist(dataset, quo_text(new_var))

  base <- dataset %>%
    filter(ABLFL == "Y") %>%
    select(!!!by_vars, !!new_var := !!source_var)

  signal_duplicate_records(base, by_vars, "Dataset contains multiple baseline records.")

  left_join(dataset, base, by = vars2chr(by_vars))
}
