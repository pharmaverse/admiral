#' Create Summary Records
#'
#' @description
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
#' @param dataset A data frame.
#'
#' @param by_vars Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [vars()] will create a
#'   groupwise summary and generate summary records for the specified groups.
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
#' @param analysis_var Analysis variable.
#'
#' @param summary_fun Function that takes as an input the `analysis_var` and
#'   performs the calculation.
#'   This can include built-in functions as well as user defined functions,
#'   for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
#'
#' @param set_values_to A list of variable name-value pairs. Use this argument
#'   if you need to change the values of any newly derived records.
#'
#'   Set a list of variables to some specified value for the new observation(s)
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a symbol, a numeric
#'   value or NA.
#'   (e.g.  `vars(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")`).
#'   More general expression are not allowed.
#'
#' @author Pavan Kumar, updated by Alana Harris
#'
#' @return A data frame of derived records.
#'
#' @family der_gen
#'
#' @keywords der_gen
#'
#' @seealso `derive_summary_records()`
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adeg <- tribble(
#'   ~USUBJID, ~EGSEQ, ~PARAM, ~AVISIT, ~EGDTC, ~AVAL, ~TRTA,
#'   "XYZ-1001", 1, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385, "",
#'   "XYZ-1001", 2, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399, "",
#'   "XYZ-1001", 3, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396, "",
#'   "XYZ-1001", 4, "QTcF Int. (msec)", "Visit 2", "2016-03-08T09:45", 384, "Placebo",
#'   "XYZ-1001", 5, "QTcF Int. (msec)", "Visit 2", "2016-03-08T09:48", 393, "Placebo",
#'   "XYZ-1001", 6, "QTcF Int. (msec)", "Visit 2", "2016-03-08T09:51", 388, "Placebo",
#'   "XYZ-1001", 7, "QTcF Int. (msec)", "Visit 3", "2016-03-22T10:45", 385, "Placebo",
#'   "XYZ-1001", 8, "QTcF Int. (msec)", "Visit 3", "2016-03-22T10:48", 394, "Placebo",
#'   "XYZ-1001", 9, "QTcF Int. (msec)", "Visit 3", "2016-03-22T10:51", 402, "Placebo",
#'   "XYZ-1002", 1, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399, "",
#'   "XYZ-1002", 2, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410, "",
#'   "XYZ-1002", 3, "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392, "",
#'   "XYZ-1002", 4, "QTcF Int. (msec)", "Visit 2", "2016-03-06T09:50", 401, "Active 20mg",
#'   "XYZ-1002", 5, "QTcF Int. (msec)", "Visit 2", "2016-03-06T09:53", 407, "Active 20mg",
#'   "XYZ-1002", 6, "QTcF Int. (msec)", "Visit 2", "2016-03-06T09:56", 400, "Active 20mg",
#'   "XYZ-1002", 7, "QTcF Int. (msec)", "Visit 3", "2016-03-24T10:50", 412, "Active 20mg",
#'   "XYZ-1002", 8, "QTcF Int. (msec)", "Visit 3", "2016-03-24T10:53", 414, "Active 20mg",
#'   "XYZ-1002", 9, "QTcF Int. (msec)", "Visit 3", "2016-03-24T10:56", 402, "Active 20mg",
#' )
#'
#' # Summarize the average of the triplicate ECG interval values (AVAL)
#' get_summary_records(
#'   adeg,
#'   by_vars = vars(USUBJID, PARAM, AVISIT),
#'   analysis_var = AVAL,
#'   summary_fun = function(x) mean(x, na.rm = TRUE),
#'   set_values_to = vars(DTYPE = "AVERAGE")
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~VSSEQ, ~PARAM, ~AVAL, ~VSSTRESU, ~VISIT, ~VSDTC,
#'   "XYZ-001-001", 1164, "Weight", 99, "kg", "Screening", "2018-03-19",
#'   "XYZ-001-001", 1165, "Weight", 101, "kg", "Run-In", "2018-03-26",
#'   "XYZ-001-001", 1166, "Weight", 100, "kg", "Baseline", "2018-04-16",
#'   "XYZ-001-001", 1167, "Weight", 94, "kg", "Week 24", "2018-09-30",
#'   "XYZ-001-001", 1168, "Weight", 92, "kg", "Week 48", "2019-03-17",
#'   "XYZ-001-001", 1169, "Weight", 95, "kg", "Week 52", "2019-04-14",
#' )
#'
#' # Set new values to any variable. Here, `DTYPE = MAXIMUM` refers to `max()` records
#' # and `DTYPE = AVERAGE` refers to `mean()` records.
#' get_summary_records(
#'   advs,
#'   by_vars = vars(USUBJID, PARAM),
#'   analysis_var = AVAL,
#'   summary_fun = max,
#'   set_values_to = vars(DTYPE = "MAXIMUM")
#' ) %>%
#'   get_summary_records(
#'     by_vars = vars(USUBJID, PARAM),
#'     analysis_var = AVAL,
#'     summary_fun = mean,
#'     set_values_to = vars(DTYPE = "AVERAGE")
#'   )
#'
#' # Sample ADEG dataset with triplicate record for only AVISIT = 'Baseline'
#' adeg <- tribble(
#'   ~USUBJID, ~EGSEQ, ~PARAM, ~AVISIT, ~EGDTC, ~AVAL, ~TRTA,
#'   "XYZ-1001", 1, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385, "",
#'   "XYZ-1001", 2, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399, "",
#'   "XYZ-1001", 3, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396, "",
#'   "XYZ-1001", 4, "QTcF Int. (msec)", "Visit 2", "2016-03-08T09:48", 393, "Placebo",
#'   "XYZ-1001", 5, "QTcF Int. (msec)", "Visit 2", "2016-03-08T09:51", 388, "Placebo",
#'   "XYZ-1001", 6, "QTcF Int. (msec)", "Visit 3", "2016-03-22T10:48", 394, "Placebo",
#'   "XYZ-1001", 7, "QTcF Int. (msec)", "Visit 3", "2016-03-22T10:51", 402, "Placebo",
#'   "XYZ-1002", 1, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399, "",
#'   "XYZ-1002", 2, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410, "",
#'   "XYZ-1002", 3, "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392, "",
#'   "XYZ-1002", 4, "QTcF Int. (msec)", "Visit 2", "2016-03-06T09:53", 407, "Active 20mg",
#'   "XYZ-1002", 5, "QTcF Int. (msec)", "Visit 2", "2016-03-06T09:56", 400, "Active 20mg",
#'   "XYZ-1002", 6, "QTcF Int. (msec)", "Visit 3", "2016-03-24T10:53", 414, "Active 20mg",
#'   "XYZ-1002", 7, "QTcF Int. (msec)", "Visit 3", "2016-03-24T10:56", 402, "Active 20mg",
#' )
#'
#' # Compute the average of AVAL only if there are more than 2 records within the
#' # by group
#' get_summary_records(
#'   adeg,
#'   by_vars = vars(USUBJID, PARAM, AVISIT),
#'   filter = n() > 2,
#'   analysis_var = AVAL,
#'   summary_fun = function(x) mean(x, na.rm = TRUE),
#'   set_values_to = vars(DTYPE = "AVERAGE")
#' )
get_summary_records <- function(dataset,
                                by_vars,
                                filter = NULL,
                                analysis_var,
                                summary_fun,
                                set_values_to = NULL) {
  assert_vars(by_vars)
  analysis_var <- assert_symbol(enquo(analysis_var))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_s3_class(summary_fun, "function")
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, analysis_var),
    check_is_grouped = FALSE
  )
  if (!is.null(set_values_to)) {
    assert_varval_list(set_values_to, optional = TRUE)
  }

  # Summarise the analysis value
  dataset %>%
    group_by(!!!by_vars) %>%
    filter_if(filter) %>%
    summarise(!!analysis_var := summary_fun(!!analysis_var)) %>%
    mutate(!!!set_values_to) %>%
    ungroup()
}
