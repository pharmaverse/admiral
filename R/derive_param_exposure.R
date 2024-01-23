#' Add an Aggregated Parameter and Derive the Associated Start and End Dates
#'
#' Add a record computed from the aggregated analysis value of another parameter and compute the
#' start (`ASTDT(M)`)and end date (`AENDT(M)`) as the minimum and maximum date by `by_vars`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified for `by_vars`, `analysis_var`, `PARAMCD`,
#'   alongside either `ASTDTM` and `AENDTM` or `ASTDT` and `AENDT` are also expected.
#'   Observations from the specified dataset are going to be used to calculate and added
#'   as new records to the input dataset (`dataset`).
#'
#'
#' @param filter
#'
#'  `r lifecycle::badge("deprecated")` Please use `filter_add` instead.
#'
#'   Filter condition as logical expression to apply during
#'   summary calculation. By default, filtering expressions are computed within
#'   `by_vars` as this will help when an aggregating, lagging, or ranking
#'   function is involved.
#'
#'   For example,
#'
#'   + `filter = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all `AVAL`
#'   values greater than mean of `AVAL` with in `by_vars`.
#'   + `filter = (dplyr::n() > 2)` will filter n count of `by_vars` greater
#'   than 2.
#'
#' @param filter_add Filter condition as logical expression to apply during
#'   summary calculation. By default, filtering expressions are computed within
#'   `by_vars` as this will help when an aggregating, lagging, or ranking
#'   function is involved.
#'
#'   For example,
#'
#'   + `filter_add = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all `AVAL`
#'   values greater than mean of `AVAL` with in `by_vars`.
#'   + `filter_add = (dplyr::n() > 2)` will filter n count of `by_vars` greater
#'   than 2.
#'
#' @param input_code Required parameter code
#'
#' The observations where `PARAMCD` equals the specified value are considered to compute the
#' summary record.
#'
#'   *Permitted Values:* A character of `PARAMCD` value
#'
#' @param analysis_var Analysis variable.
#'
#' @param summary_fun Function that takes as an input the `analysis_var` and
#'   performs the calculation.
#'   This can include built-in functions as well as user defined functions,
#'   for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param set_values_to Variable-value pairs
#'
#'   Set a list of variables to some specified value for the new observation(s)
#'   + LHS refer to a variable. It is expected that at least `PARAMCD` is defined.
#'   + RHS refers to the values to set to the variable. This can be a string, a symbol, a numeric
#'   value, `NA`, or an expression.
#'   (e.g.  `exprs(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")`).
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @details For each group (with respect to the variables specified for the `by_vars` parameter),
#' an observation is added to the output dataset and the defined values are set to the defined
#' variables
#'
#' @return The input dataset with a new record added for each group (with respect to the variables
#' specified for the `by_vars` parameter). That is, a variable will only
#' be populated in this new record if it is specified in `by_vars`.
#' For each new record,
#' + the variable specified `analysis_var` is computed as defined by `summary_fun`,
#' + the variable(s) specified on the LHS of `set_values_to` are set to their paired value (RHS).
#' In addition, the start and end date are computed as the minimum/maximum dates by `by_vars`.
#'
#' If the input datasets contains
#' + both `AxxDTM` and `AxxDT` then all `ASTDTM`,`AENDTM`, `ASTDT`, `AENDT` are computed
#' + only `AxxDTM` then `ASTDTM`,`AENDTM` are computed
#' + only `AxxDT` then `ASTDT`,`AENDT` are computed.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#' library(stringr, warn.conflicts = FALSE)
#' adex <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC, ~VISIT, ~ASTDT, ~AENDT,
#'   "1015", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "1015", "DOSE", 85, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "1015", "DOSE", 82, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "1015", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "1015", "ADJ", NA, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "1015", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "1017", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "1017", "DOSE", 50, NA_character_, "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "1017", "DOSE", 65, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02"),
#'   "1017", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "1017", "ADJ", NA, "ADVERSE EVENT", "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "1017", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02")
#' ) %>%
#'   mutate(ASTDTM = ymd_hms(paste(ASTDT, "00:00:00")), AENDTM = ymd_hms(paste(AENDT, "00:00:00")))
#'
#' # Cumulative dose
#' adex %>%
#'   derive_param_exposure(
#'     dataset_add = adex,
#'     by_vars = exprs(USUBJID),
#'     set_values_to = exprs(
#'       PARAMCD = "TDOSE",
#'       PARCAT1 = "OVERALL",
#'       AVAL = sum(AVAL, na.rm = TRUE)
#'     ),
#'     input_code = "DOSE"
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
#'
#' # average dose in w2-24
#' adex %>%
#'   derive_param_exposure(
#'     dataset_add = adex,
#'     by_vars = exprs(USUBJID),
#'     filter_add = VISIT %in% c("WEEK 2", "WEEK 24"),
#'     set_values_to = exprs(
#'       PARAMCD = "AVDW224",
#'       PARCAT1 = "WEEK2-24",
#'       AVAL = mean(AVAL, na.rm = TRUE)
#'     ),
#'     input_code = "DOSE"
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
#'
#' # Any dose adjustment?
#' adex %>%
#'   derive_param_exposure(
#'     dataset_add = adex,
#'     by_vars = exprs(USUBJID),
#'     set_values_to = exprs(
#'       PARAMCD = "TADJ",
#'       PARCAT1 = "OVERALL",
#'       AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
#'     ),
#'     input_code = "ADJ"
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
derive_param_exposure <- function(dataset = NULL,
                                  dataset_add,
                                  by_vars,
                                  input_code,
                                  analysis_var,
                                  summary_fun,
                                  filter = NULL,
                                  filter_add = NULL,
                                  set_values_to = NULL) {
  by_vars <- assert_vars(by_vars)
  if (!missing(analysis_var) || !missing(summary_fun)) {
    deprecate_stop(
      "1.1.0",
      I("derive_param_exposure(anaylsis_var = , summary_fun = )"),
      "derive_param_exposure(set_values_to = )"
    )
  }

  dtm <- c("ASTDTM", "AENDTM") %in% colnames(dataset)
  dt <- c("ASTDT", "AENDT") %in% colnames(dataset)
  set_dtm <- NULL
  set_dt <- NULL
  if (all(dtm)) {
    dates <- exprs(ASTDTM, AENDTM)
    set_dtm <- exprs(
      ASTDTM = min(ASTDTM, na.rm = TRUE),
      AENDTM = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
    )
  } else {
    dates <- exprs(ASTDT, AENDT)
  }
  if (all(dt)) {
    set_dt <- exprs(
      ASTDT = min(ASTDT, na.rm = TRUE),
      AENDT = max(coalesce(AENDT, ASTDT), na.rm = TRUE)
    )
  }

  assert_data_frame(dataset, required_vars = by_vars, optional = TRUE)
  assert_data_frame(dataset_add,
    required_vars = expr_c(by_vars, exprs(PARAMCD), dates)
  )

  if (!missing(filter)) {
    deprecate_stop(
      "1.1.0",
      I("derive_param_exposure(filter = )"),
      "derive_param_exposure(filter_add = )"
    )
    filter_add <- assert_filter_cond(enexpr(filter), optional = TRUE)
  }
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(input_code)
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(input_code, values = params_available)

  if (is.null(filter_add)) {
    filter_add <- TRUE
  }

  derive_summary_records(
    dataset,
    dataset_add,
    by_vars = by_vars,
    filter_add = PARAMCD == !!input_code & !!filter_add,
    set_values_to = exprs(
      !!!set_dtm,
      !!!set_dt,
      !!!set_values_to
    )
  )
}
