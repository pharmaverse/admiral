#' Create Summary Records
#'
#' @description
#'
#' `r lifecycle::badge("deprecated")` The `get_summary_records()` has been
#' deprecated in favor of `derive_summary_records()` (call it with the `dataset_add`
#' argument and without the `dataset` argument).
#'
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple records. The ADaM basic dataset
#' structure variable `DTYPE` is available to indicate when a new derived
#' records has been added to a dataset.
#'
#' @details
#' This function only creates derived observations and does not append them
#' to the original dataset observations. If you would like to this instead,
#' see the `derive_summary_records()` function.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param by_vars Grouping variables
#'
#'   Variables to consider for generation of groupwise summary records.
#'
#'  `r roxygen_param_by_vars()`
#'
#' @param filter Filter condition as logical expression to apply during
#'   summary calculation. By default, filtering expressions are computed within
#'   `by_vars` as this will help when an aggregating, lagging, or ranking
#'   function is involved.
#'
#'   For example,
#'
#'   + `filter_rows = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all AVAL
#'   values greater than mean of AVAL with in `by_vars`.
#'   + `filter_rows = (dplyr::n() > 2)` will filter n count of `by_vars` greater
#'   than 2.
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations.
#'
#'   Set a list of variables to some specified value for the new records
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value, an expression or NA. If summary functions are
#'   used, the values are summarized by the variables specified for `by_vars`.
#'
#'   For example:
#'   ```
#'     set_values_to = exprs(
#'       AVAL = sum(AVAL),
#'       PARAMCD = "TDOSE",
#'       PARCAT1 = "OVERALL"
#'     )
#'   ```
#'
#' @return A data frame of derived records.
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [derive_summary_records()], [derive_var_merged_summary()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adeg <- tribble(
#'   ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
#'   "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
#'   "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
#'   "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
#'   "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:45", 384,   "Placebo",
#'   "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
#'   "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
#'   "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:45", 385,   "Placebo",
#'   "XYZ-1001", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
#'   "XYZ-1001", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
#'   "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
#'   "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
#'   "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
#'   "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:50", 401,   "Active 20mg",
#'   "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
#'   "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
#'   "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:50", 412,   "Active 20mg",
#'   "XYZ-1002", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
#'   "XYZ-1002", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
#' )
#'
#' # Summarize the average of the triplicate ECG interval values (AVAL)
#' get_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' )
#'
#' # Derive more than one summary variable
#' get_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL),
#'     ASTDTM = min(convert_dtc_to_dtm(EGDTC)),
#'     AENDTM = max(convert_dtc_to_dtm(EGDTC)),
#'     DTYPE = "AVERAGE"
#'   )
#' )
#'
#' # Sample ADEG dataset with triplicate record for only AVISIT = 'Baseline'
#' adeg <- tribble(
#'   ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
#'   "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
#'   "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
#'   "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
#'   "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
#'   "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
#'   "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
#'   "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
#'   "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
#'   "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
#'   "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
#'   "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
#'   "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
#'   "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
#'   "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
#' )
#'
#' # Compute the average of AVAL only if there are more than 2 records within the
#' # by group
#' get_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   filter = n() > 2,
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' )
get_summary_records <- function(dataset,
                                by_vars,
                                filter = NULL,
                                set_values_to = NULL) {
  deprecate_inform(
    when = "1.2.0",
    what = "get_summary_records()",
    with = "derive_summary_records()",
    details = c(
      x = "This message will turn into a warning at the beginning of 2026.",
      i = "See admiral's deprecation guidance:
      https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )

  assert_vars(by_vars)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = by_vars,
    check_is_grouped = FALSE
  )
  assert_varval_list(set_values_to)

  # Summarise the analysis value
  dataset %>%
    group_by(!!!by_vars) %>%
    filter_if(filter) %>%
    summarise(!!!set_values_to) %>%
    ungroup()
}
