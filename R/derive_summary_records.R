#' Add New Records Within By Groups Using Aggregation Functions
#'
#' @description
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple records. The ADaM basic dataset
#' structure variable `DTYPE` is available to indicate when a new derived
#' records has been added to a dataset, if the derivation deviates from the
#' standard derivation of the parameter.
#'
#' @details
#' For the newly derived records, only variables specified within `by_vars` or
#' `set_values_to` will be populated. All other variables will be set to `NA`.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#'   If the argument is not specified (or set to `NULL`), a new dataset is
#'   created. Otherwise, the new records are appended to the specified dataset.
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified for `by_vars` are expected.
#'   Observations from the specified dataset are going to be used to calculate and added
#'   as new records to the input dataset (`dataset`).
#'
#' @permitted [dataset]
#'
#' @param dataset_ref Reference dataset
#'
#'   The variables specified for `by_vars` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [exprs()] will create a
#'   groupwise summary and generate summary records for the specified groups.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
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
#' @permitted [condition]
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
#' @permitted [expr_list_summary]
#'
#' @param missing_values Values for missing summary values
#'
#'   For observations of the reference dataset (`dataset_ref`) which do not have a
#'   complete mapping defined by the summarization defined in `set_values_to`.  Only variables
#'   specified for `set_values_to` can be specified for `missing_values`.
#'
#' @permitted [expr_list_summary]
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @seealso [derive_var_merged_summary()]
#'
#' @export
#'
#' @examplesx
#'
#' @caption Data setup
#'
#' @info The following examples use the ECG dataset below as a basis.
#'
#' @code
#' library(tibble, warn.conflicts = FALSE)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adeg <- tribble(
#'   ~USUBJID,   ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL,
#'   "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,
#'   "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,
#'   "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,
#'   "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,
#'   "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 200,
#'   "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,
#'   "XYZ-1002", "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,
#'   "XYZ-1002", "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402
#' ) %>%
#'   mutate(ADTM = convert_dtc_to_dtm(EGDTC))
#'
#' @caption Summarize one or more variables using summary functions
#'
#' @info A derived record is generated for each subject, containing the mean of the triplicate ECG
#' interval values (`AVAL`) and the latest measurement's time (`ADTM`) by using summary functions
#' within the `set_values_to` argument.
#'
#' @code
#' derive_summary_records(
#'   adeg,
#'   dataset_add = adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     ADTM = max(ADTM),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
#'
#' @info Functions such as `all()` and `any()` are also often useful when creating
#' summary records. For instance, the above example can be extended to flag which derived
#' records were affected by outliers. Note that the outlier flag is created before `AVAL`
#' is set for the summary record. Otherwise, referencing `AVAL` later on would pick up the
#' `AVAL` from the summary record rather than the source records.
#'
#' @code
#' derive_summary_records(
#'   adeg,
#'   dataset_add = adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     OUTLIEFL = if_else(any(AVAL >= 500 | AVAL <= 300), "Y", "N"),
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     ADTM = max(ADTM),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
#'
#' @caption Restricting source records (`filter_add`)
#'
#' @info The `filter_add` argument can be used to restrict the records that are being summarized.
#' For instance, the mean of the triplicates above can be computed only for the baseline records
#' by passing `filter_add = AVISIT == "Baseline"`.
#'
#' @code
#' derive_summary_records(
#'   adeg,
#'   dataset_add = adeg,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   filter_add = AVISIT == "Baseline",
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     DTYPE = "AVERAGE"
#'   )
#' ) %>%
#'   arrange(USUBJID, AVISIT)
#'
#' @info Summary functions can also be used within `filter_add` to filter based on conditions
#' applied to the whole of the by group specified in `by_vars`. For instance, the mean of
#' the triplicates can be computed only for by groups  which do indeed contain three records
#' by passing `filter_add = n() > 2`.
#'
#' @code
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
#'
#' @caption Adding records for groups not in source (`dataset_ref` and `missing_values`)
#'
#' @info Adding records for groups which are not in the source data can be achieved by
#' specifying a reference dataset in the `dataset_ref` argument. For example, specifying
#' the input dataset `adeg_allparamvis` (containing an extra `"Visit 2"` for patient
#' `1002`) ensures a summary record is derived for that visit as well. For these records,
#' the values of the analysis variables to be populated should be specified within the
#' `missing_values` argument. Here, `DTYPE = "PHANTOM"` was chosen as `AVAL` is set to
#' missing.
#'
#' @code
#' adeg_allparamvis <- tribble(
#'   ~USUBJID,   ~PARAM,             ~AVISIT,
#'   "XYZ-1001", "QTcF Int. (msec)", "Baseline",
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 2",
#'   "XYZ-1001", "QTcF Int. (msec)", "Visit 3",
#'   "XYZ-1002", "QTcF Int. (msec)", "Baseline",
#'   "XYZ-1002", "QTcF Int. (msec)", "Visit 2",
#'   "XYZ-1002", "QTcF Int. (msec)", "Visit 3"
#' )
#'
#' derive_summary_records(
#'   adeg,
#'   dataset_add = adeg,
#'   dataset_ref = adeg_allparamvis,
#'   by_vars = exprs(USUBJID, PARAM, AVISIT),
#'   set_values_to = exprs(
#'     AVAL = mean(AVAL, na.rm = TRUE),
#'     ADTM = max(ADTM),
#'     DTYPE = "AVERAGE"
#'   ),
#'   missing_values = exprs(
#'     AVAL = NA,
#'     ADTM = NA,
#'     DTYPE = "PHANTOM"
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
  assert_data_frame(dataset, optional = TRUE)
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
    new_ref_obs <- anti_join(
      select(dataset_ref, !!!by_vars),
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
