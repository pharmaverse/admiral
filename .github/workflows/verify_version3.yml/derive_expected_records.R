#' Derive Expected Records
#'
#' Add expected records as new observations for each 'by group' when the dataset
#' contains missing observations.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("dataset_ref", "by_vars"))`
#'
#' @param dataset_ref Expected observations dataset
#'
#'   Data frame with the expected observations, e.g., all the expected
#'   combinations of `PARAMCD`, `PARAM`, `AVISIT`, `AVISITN`, ...
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` those observations from `dataset_ref`
#'   are added to the output dataset which do not have a corresponding observation
#'   in the input dataset.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations.
#'
#'   A list of variable name-value pairs is expected.
#'   + LHS refers to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value, `NA`, or expressions, e.g., `exprs(PARAMCD =
#'   "TDOSE", PARCAT1 = "OVERALL")`.
#'
#' @details For each group (the variables specified in the `by_vars` parameter),
#' those records from `dataset_ref` that are missing in the input
#' dataset are added to the output dataset.
#'
#' @return The input dataset with the missed expected observations added for each
#' `by_vars`. Note, a variable will only be populated in the new parameter rows
#' if it is specified in `by_vars` or `set_values_to`.
#'
#' @keywords der_prm_bds_findings
#' @family der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adqs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
#'   "1",      "a",             1, "WEEK 1",   10,
#'   "1",      "b",             1, "WEEK 1",   11,
#'   "2",      "a",             2, "WEEK 2",   12,
#'   "2",      "b",             2, "WEEK 2",   14
#' )
#'
#' # Example 1. visit variables are parameter independent
#' parm_visit_ref <- tribble(
#'   ~AVISITN, ~AVISIT,
#'   1,        "WEEK 1",
#'   2,        "WEEK 2"
#' )
#'
#' derive_expected_records(
#'   dataset = adqs,
#'   dataset_ref = parm_visit_ref,
#'   by_vars = exprs(USUBJID, PARAMCD),
#'   set_values_to = exprs(DTYPE = "DERIVED")
#' )
#'
#' # Example 2. visit variables are parameter dependent
#' parm_visit_ref <- tribble(
#'   ~PARAMCD, ~AVISITN, ~AVISIT,
#'   "a",             1, "WEEK 1",
#'   "a",             2, "WEEK 2",
#'   "b",             1, "WEEK 1"
#' )
#'
#' derive_expected_records(
#'   dataset = adqs,
#'   dataset_ref = parm_visit_ref,
#'   by_vars = exprs(USUBJID, PARAMCD),
#'   set_values_to = exprs(DTYPE = "DERIVED")
#' )
#'
derive_expected_records <- function(dataset,
                                    dataset_ref,
                                    by_vars = NULL,
                                    set_values_to = NULL) {
  # Check input parameters
  assert_vars(by_vars, optional = TRUE)
  assert_data_frame(dataset_ref)
  assert_data_frame(
    dataset,
    required_vars = expr_c(by_vars, chr2vars(colnames(dataset_ref)))
  )
  assert_varval_list(set_values_to, optional = TRUE)

  # Derive expected records
  ## ids: Variables from by_vars but not in dataset_ref
  ids <- dataset %>%
    select(!!!setdiff(by_vars, chr2vars(colnames(dataset_ref)))) %>%
    distinct()

  if (ncol(ids) > 0) {
    exp_obsv <- ids %>%
      crossing(dataset_ref)
  } else {
    exp_obsv <- dataset_ref
  } # tmp workaround, update after using tidyr 1.2.0

  exp_obs_vars <- exp_obsv %>%
    colnames()

  new_obs <- exp_obsv %>%
    mutate(!!!set_values_to) %>%
    anti_join(dataset %>% select(all_of(exp_obs_vars)), by = c(exp_obs_vars))

  # Combine dataset + newly added records
  bind_rows(dataset, new_obs) %>%
    arrange(!!!chr2vars(exp_obs_vars))
}
