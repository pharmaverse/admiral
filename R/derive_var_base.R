#' Derive Baseline Variables
#'
#' Derive baseline variables, e.g. `BASE` or `BNRIND`, in a BDS dataset
#'
#' @param dataset The input dataset
#'
#' @param by_vars Grouping variables uniquely identifying a set
#'        of records for which to calculate `new_var`
#'
#' @param source_var The column from which to extract the baseline value, e.g. `AVAL`
#'
#' @param new_var The name of the newly created baseline column, e.g. `BASE`
#'
#' @param filter The condition used to filter `dataset` for baseline records.
#'
#'   By default `ABLFL == "Y"`
#'
#' @return
#' A new `data.frame` containing all records and variables of the input
#' dataset plus the `new_var` variable
#'
#' @details
#' For each `by_vars` group the baseline record is identified by filtering using the
#' condition specified by `filter` which defaults to `ABLFL == "Y"`. Subsequently,
#' every value of the `new_var` variable for the `by_vars` group is set to the
#' value of the `source_var` variable of the baseline record. In case there are
#' multiple baseline records within `by_vars` an error is issued.
#'
#' @export
#'
#' @author Thomas Neitmann
#'
#' @keywords bds derivation
#'
#' @examples
#' dataset <- tibble::tribble(
#'   ~STUDYID, ~USUBJID,   ~PARAMCD, ~AVAL,   ~AVALC,      ~AVISIT, ~ABLFL,
#'   "TEST01",  "PAT01",  "PARAM01", 10.12,       NA,   "Baseline",    "Y",
#'   "TEST01",  "PAT01",  "PARAM01", 9.700,       NA,      "Day 7",    "N",
#'   "TEST01",  "PAT01",  "PARAM01", 15.01,       NA,     "Day 14",    "N",
#'   "TEST01",  "PAT01",  "PARAM02", 8.350,       NA,   "Baseline",    "Y",
#'   "TEST01",  "PAT01",  "PARAM02",    NA,       NA,      "Day 7",    "N",
#'   "TEST01",  "PAT01",  "PARAM02", 8.350,       NA,     "Day 14",    "N",
#'   "TEST01",  "PAT01",  "PARAM03",    NA,    "LOW",   "Baseline",    "Y",
#'   "TEST01",  "PAT01",  "PARAM03",    NA,    "LOW",      "Day 7",    "N",
#'   "TEST01",  "PAT01",  "PARAM03",    NA, "MEDIUM",     "Day 14",    "N",
#'   "TEST01",  "PAT01",  "PARAM04",    NA,   "HIGH",   "Baseline",    "Y",
#'   "TEST01",  "PAT01",  "PARAM04",    NA,   "HIGH",      "Day 7",    "N",
#'   "TEST01",  "PAT01",  "PARAM04",    NA, "MEDIUM",     "Day 14",    "N"
#' )
#'
#' ## Derive `BASE` variable from `AVAL`
#' derive_var_base(
#'   dataset,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   source_var = AVAL,
#'   new_var = BASE
#' )
#'
#' ## Derive `BASEC` variable from `AVALC`
#' derive_var_base(
#'   dataset,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   source_var = AVALC,
#'   new_var = BASEC
#' )
#'
#' ## Derive `BNRIND` variable from `ANRIND`
#' if (FALSE) {
#'   derive_var_base(
#'     dataset,
#'     by_vars = vars(USUBJID, PARAMCD),
#'     source_var = ANRIND,
#'     new_var = BNRIND
#'   )
#' }
derive_var_base <- function(dataset,
                            by_vars,
                            source_var = AVAL,
                            new_var = BASE,
                            filter = ABLFL == "Y") {
  by_vars <- assert_vars(by_vars)
  source_var <- assert_symbol(enquo(source_var))
  new_var <- assert_symbol(enquo(new_var))
  filter <- assert_filter_cond(enquo(filter))
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, source_var)
  )
  warn_if_vars_exist(dataset, quo_text(new_var))

  derive_vars_merged(
    dataset,
    dataset_add = dataset,
    filter_add = !!filter,
    new_vars = vars(!!new_var := !!source_var),
    by_vars = by_vars,
    duplicate_msg = paste(
      "Input dataset contains multiple baseline records with respect to",
      enumerate(vars2chr(by_vars))
    )
  )
}
