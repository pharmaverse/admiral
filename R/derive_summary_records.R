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
#' @param dataset  `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_add Additional dataset
#'
#'    The variables specified for `by_vars` are expected.
#'   Observations from the specified dataset are going to be used to calculate and added
#'   as new records to the input dataset (`dataset`).
#'
#' @param dataset_ref Reference dataset
#'
#'   The variables specified for `by_vars` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param by_vars Grouping variables
#'
#'   Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [exprs()] will create a
#'   groupwise summary and generate summary records for the specified groups.
#'
#'   `r roxygen_param_by_vars()`
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
#'       DTYPE = "AVERAGE",
#'     )
#'   ```
#'
#' @param missing_values Values for missing summary values
#'
#'   For observations of the reference dataset (`dataset_ref`) which do not have a
#'   complete mapping defined by the summarization defined in `set_values_to`.  Only variables
#'   specified for `set_values_to` can be specified for `missing_values`.
#'
#'   *Permitted Values*: named list of expressions, e.g.,
#'   `exprs(AVAL = -9999)`
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
#'   dataset_add = adeg,
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
#'   dataset_add = adeg,
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
#'   dataset_add = adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   filter_add = n() > 2,
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
derive_summary_records <- function(dataset = NULL,
                                   dataset_add,
                                   dataset_ref = NULL,
                                   by_vars,
                                   filter_add = NULL,
                                   set_values_to,
                                   missing_values = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = by_vars, optional = TRUE)
  assert_data_frame(dataset_add, required_vars = by_vars)
  assert_data_frame(
    dataset_ref,
    required_vars = by_vars,
    optional = TRUE
  )

  assert_varval_list(set_values_to)
  assert_expr_list(missing_values, named = TRUE, optional = TRUE)

  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)

  summary_records <- dataset_add %>%
    group_by(!!!by_vars) %>%
    filter_if(filter_add) %>%
    summarise(!!!set_values_to) %>%
    ungroup()

  df_return <- bind_rows(
    dataset,
    summary_records
  )

  if (!is.null(dataset_ref)) {
    add_vars <- colnames(dataset_add)
    ref_vars <- colnames(dataset_ref)

    new_ref_obs <- anti_join(
      select(dataset_ref, intersect(add_vars, ref_vars)),
      select(summary_records, !!!by_vars),
      by = map_chr(by_vars, as_name)
    )

    if (!is.null(missing_values)) {
      new_ref_obs <- new_ref_obs %>%
        mutate(!!!missing_values)
    }

    df_return <- bind_rows(
      df_return,
      new_ref_obs
    )
  }

  df_return
}
