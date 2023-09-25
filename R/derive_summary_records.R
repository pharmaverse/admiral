#' Add New Records Within By Groups Using Aggregation Functions
#'
#' @description
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple records. The ADaM basic dataset
#' structure variable `DTYPE` is available to indicate when a new derived
#' records has been added to a dataset.
#'
#' @details
#' When all records have same values within `by_vars` then this function will
#' retain those common values in the newly derived records. Otherwise new value
#' will be set to `NA`.
#'
#' @param dataset A data frame.
#'
#' @param by_vars Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [exprs()] will create a
#'   groupwise summary and generate summary records for the specified groups.
#'
#' @param filter Filter condition as logical expression to apply during
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
#' @inheritParams get_summary_records
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @seealso [get_summary_records()], [derive_var_merged_summary()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
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
#' ) %>%
#'   mutate(
#'     ADTM = convert_dtc_to_dtm(EGDTC)
#'   )
#'
#' # Summarize the average of the triplicate ECG interval values (AVAL)
#' derive_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
#'
#' # Derive more than one summary variable
#' derive_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL),
#'     ADTM = max(ADTM),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT) %>%
#'   select(-EGSEQ, -TRTA)
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
#' derive_summary_records(
#'   adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   filter = n() > 2,
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
derive_summary_records <- function(dataset,
                                   by_vars,
                                   filter = NULL,
                                   analysis_var,
                                   summary_fun,
                                   set_values_to) {
  assert_vars(by_vars)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = by_vars
  )
  assert_varval_list(set_values_to)

  if (!missing(analysis_var) || !missing(summary_fun)) {
    deprecate_warn(
      "1.0.0",
      I("derive_summary_records(anaylsis_var = , summary_fun = )"),
      "derive_summary_records(set_values_to = )"
    )
    analysis_var <- assert_symbol(enexpr(analysis_var))
    assert_s3_class(summary_fun, "function")
    set_values_to <- exprs(!!analysis_var := {{ summary_fun }}(!!analysis_var), !!!set_values_to)
  }

  # Summarise the analysis value and bind to the original dataset
  bind_rows(
    dataset,
    get_summary_records(
      dataset,
      by_vars = by_vars,
      filter = !!filter,
      set_values_to = set_values_to
    )
  )
}
