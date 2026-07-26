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
#' @seealso [derive_summary_records()], [derive_vars_merged_summary()]
#'
#' @export
#'
get_summary_records <- function(dataset,
                                by_vars,
                                filter = NULL,
                                set_values_to = NULL) {
  deprecate_stop(
    when = "1.6.0",
    what = "get_summary_records()",
    with = "derive_summary_records()"
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
