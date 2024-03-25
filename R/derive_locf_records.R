#' Derive LOCF (Last Observation Carried Forward) Records
#'
#' Adds LOCF records as new observations for each 'by group' when the dataset
#' does not contain observations for missed visits/time points.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars", "analysis_var", "order", "keep_vars"))`
#'
#' @param dataset_ref Expected observations dataset
#'
#'   Data frame with all the combinations of `PARAMCD`, `PARAM`, `AVISIT`,
#'   `AVISITN`, ... which are expected in the dataset is expected.
#'
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` those observations from `dataset_ref`
#'   are added to the output dataset which do not have a corresponding observation
#'   in the input dataset or for which `analysis_var` is `NA` for the corresponding observation
#'   in the input dataset.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param analysis_var Analysis variable.
#'
#'   *Default*: `AVAL`
#'
#'   *Permitted Values*: a variable
#'
#' @param order Sort order
#'
#'   The dataset is sorted by `order` before carrying the last observation
#'   forward (e.g. `AVAL`) within each `by_vars`.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @param keep_vars Variables that need carrying the last observation forward
#'
#'   Keep variables that need carrying the last observation forward other than `analysis_var`
#'   (e.g., `PARAMN`, `VISITNUM`). If by default `NULL`, only variables specified in
#'   `by_vars` and `analysis_var` will be populated in the newly created records.
#'
#' @author G Gayatri
#'
#' @details For each group (with respect to the variables specified for the
#' by_vars parameter) those observations from `dataset_ref` are added to
#' the output dataset
#' - which do not have a corresponding observation in the input dataset or
#' - for which `analysis_var` is NA for the corresponding observation in the input dataset.
#'
#'   For the new observations, `analysis_var` is set to the non-missing `analysis_var` of the
#'   previous observation in the input dataset (when sorted by `order`) and
#'   `DTYPE` is set to "LOCF".
#'
#' @return The input dataset with the new "LOCF" observations added for each
#' `by_vars`. Note, a variable will only be populated in the new parameter rows
#' if it is specified in `by_vars`.
#'
#' @keywords der_prm_bds_findings
#' @family der_prm_bds_findings
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(tibble)
#'
#' advs <- tribble(
#'   ~STUDYID,  ~USUBJID,      ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT,
#'   "CDISC01", "01-701-1015", "PULSE",        1,    61,        0, "BASELINE",
#'   "CDISC01", "01-701-1015", "PULSE",        1,    60,        2, "WEEK 2",
#'   "CDISC01", "01-701-1015", "DIABP",        2,    51,        0, "BASELINE",
#'   "CDISC01", "01-701-1015", "DIABP",        2,    50,        2, "WEEK 2",
#'   "CDISC01", "01-701-1015", "DIABP",        2,    51,        4, "WEEK 4",
#'   "CDISC01", "01-701-1015", "DIABP",        2,    50,        6, "WEEK 6",
#'   "CDISC01", "01-701-1015", "SYSBP",        3,   121,        0, "BASELINE",
#'   "CDISC01", "01-701-1015", "SYSBP",        3,   121,        2, "WEEK 2",
#'   "CDISC01", "01-701-1015", "SYSBP",        3,   121,        4, "WEEK 4",
#'   "CDISC01", "01-701-1015", "SYSBP",        3,   121,        6, "WEEK 6",
#'   "CDISC01", "01-701-1028", "PULSE",        1,    65,        0, "BASELINE",
#'   "CDISC01", "01-701-1028", "DIABP",        2,    79,        0, "BASELINE",
#'   "CDISC01", "01-701-1028", "DIABP",        2,    80,        2, "WEEK 2",
#'   "CDISC01", "01-701-1028", "DIABP",        2,    NA,        4, "WEEK 4",
#'   "CDISC01", "01-701-1028", "DIABP",        2,    NA,        6, "WEEK 6",
#'   "CDISC01", "01-701-1028", "SYSBP",        3,   130,        0, "BASELINE",
#'   "CDISC01", "01-701-1028", "SYSBP",        3,   132,        2, "WEEK 2"
#' )
#'
#'
#' # A dataset with all the combinations of PARAMCD, PARAM, AVISIT, AVISITN, ... which are expected.
#' advs_expected_obsv <- tribble(
#'   ~PARAMCD, ~AVISITN, ~AVISIT,
#'   "PULSE",         0, "BASELINE",
#'   "PULSE",         6, "WEEK 6",
#'   "DIABP",         0, "BASELINE",
#'   "DIABP",         2, "WEEK 2",
#'   "DIABP",         4, "WEEK 4",
#'   "DIABP",         6, "WEEK 6",
#'   "SYSBP",         0, "BASELINE",
#'   "SYSBP",         2, "WEEK 2",
#'   "SYSBP",         4, "WEEK 4",
#'   "SYSBP",         6, "WEEK 6"
#' )
#'
#' derive_locf_records(
#'   dataset = advs,
#'   dataset_ref = advs_expected_obsv,
#'   by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'   order = exprs(AVISITN, AVISIT),
#'   keep_vars = exprs(PARAMN)
#' )
#'
derive_locf_records <- function(dataset,
                                dataset_ref,
                                by_vars,
                                analysis_var = AVAL,
                                order,
                                keep_vars = NULL) {
  #### Input Checking ####
  analysis_var <- assert_symbol(enexpr(analysis_var))

  # Check if input parameters is a valid list of variables
  assert_vars(by_vars, optional = TRUE)
  assert_vars(keep_vars, optional = TRUE)
  assert_expr_list(order)

  # Check by_vars and order variables in input datasets
  assert_data_frame(dataset_ref)
  assert_data_frame(
    dataset,
    required_vars = expr_c(
      by_vars, analysis_var, extract_vars(order), keep_vars,
      chr2vars(colnames(dataset_ref))
    )
  )


  #### Prepping 'dataset_ref' ####


  # Get the IDs from input dataset for which the expected observations are to be added

  ids <- dataset %>%
    select(!!!setdiff(by_vars, chr2vars(colnames(dataset_ref)))) %>%
    distinct()

  exp_obsv <- ids %>%
    crossing(dataset_ref)



  ##### Add LOCF records ####

  # Split input dataset into the missing and non-missing analysis_var (e.g., AVAL) records
  aval_missing <- dataset %>%
    filter(is.na(!!analysis_var))

  aval_not_missing <- dataset %>%
    drop_na(!!analysis_var)


  # Get the variable names from the expected observation dataset
  exp_obs_vars <- exp_obsv %>%
    colnames()


  # Get unique combination of visits/timepoints per parameter per subject
  # from the original input dataset (with non-missing analysis_var)
  advs_unique_original <- aval_not_missing %>%
    select(all_of(exp_obs_vars)) %>%
    distinct()


  tmp_dtype <- get_new_tmp_var(exp_obsv, prefix = "tmp_dtype")

  # Get all the expected observations that are to be added to the input
  # dataset (with non-missing analysis_var)
  advs_exp_obsv3 <- exp_obsv %>%
    mutate(!!tmp_dtype := "LOCF") %>%
    anti_join(advs_unique_original, by = c(exp_obs_vars))

  # Merge the expected observations with the input dataset (with non-missing analysis_var)
  # Arrange the dataset by 'order' and group it by 'by_vars'
  # Use fill() to fill the analysis_var from the previous observation for the newly added records


  aval_not_missing_locf <- aval_not_missing %>%
    full_join(advs_exp_obsv3, by = c(exp_obs_vars))

  if ("DTYPE" %in% colnames(aval_not_missing)) {
    aval_not_missing_locf <- aval_not_missing_locf %>%
      mutate(DTYPE = if_else(!!tmp_dtype == "LOCF", "LOCF", DTYPE, missing = DTYPE)) %>%
      select(-!!tmp_dtype)
  } else {
    aval_not_missing_locf <- rename(aval_not_missing_locf, DTYPE = !!tmp_dtype)
  }

  aval_not_missing_locf <- aval_not_missing_locf %>%
    arrange(!!!by_vars, !!!order) %>%
    group_by(!!!by_vars) %>%
    fill(!!analysis_var, !!!keep_vars) %>%
    ungroup()



  # Output dataset - merge the analysis_var missing with non-missing+newly added LOCF records
  bind_rows(aval_not_missing_locf, aval_missing)
}
