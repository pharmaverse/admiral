#' Derive LOCF (Last Observation Carried Forward) Records
#'
#' Adds LOCF records as new observations for each 'by group' when the dataset
#' does not contain observations for missed visits/time points and when analysis
#' value is missing.
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
#'
#' @param id_vars_ref Grouping variables in expected observations dataset
#'
#'  The variables to group by in `dataset_ref` when determining which observations should be
#'  added to the input dataset.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @default All the variables in `dataset_ref`
#'
#'
#' @param analysis_var Analysis variable.
#'
#' @permitted a variable
#'
#' @param imputation Select the mode of imputation:
#'
#'   `add`: Keep all original records and add imputed records for missing
#'   timepoints and missing `analysis_var` values from `dataset_ref`.
#'
#'   `update`: Update records with missing `analysis_var` and add imputed records
#'   for missing timepoints from `dataset_ref`.
#'
#'   `update_add`: Keep all original records, update records with missing `analysis_var`
#'    and add imputed records for missing timepoints from `dataset_ref`.
#'
#'
#' @permitted One of these 3 values: `"add"`, `"update"`, `"update_add"`
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
#' - for which `analysis_var` is `NA` for the corresponding observation in the input dataset.
#'
#'   For the new observations, `analysis_var` is set to the non-missing `analysis_var` of the
#'   previous observation in the input dataset (when sorted by `order`) and
#'   `DTYPE` is set to "LOCF".
#'
#'   The `imputation` argument decides whether to update the existing observation when
#'   `analysis_var` is `NA` (`"update"` and `"update_add"`), or to add a new observation from
#'   `dataset_ref` instead (`"add"`).
#'
#' @return The input dataset with the new "LOCF" observations added for each
#' `by_vars`, based on the value passed to the `imputation` argument.
#'
#' @keywords der_prm_bds_findings
#' @family der_prm_bds_findings
#'
#' @export
#'
#' @examplesx
#' @caption Add records for missing analysis variable using reference dataset
#' @info Imputed records should be added for missing timepoints and for missing
#'       `analysis_var` (from `dataset_ref`), while retaining all original records.
#'
#' - The reference dataset for the imputed records is specified by the `dataset_add`
#'   argument. It should contain all expected combinations of variables. In this case,
#'   `advs_expected_obsv` is created by `crossing()` datasets `paramcd` and `avisit`, which
#'   includes all combinations of PARAMCD, AVISITN, and AVISIT.
#' - The groups for which new records are added are specified by the `by_vars`
#'   argument. Here, one record should be added for each *subject* and *parameter*.
#'   Therefore, `by_vars = exprs(STUDYID, USUBJID, PARAMCD)` is specified.
#' - The imputation method is specified using the `imputation` argument. In this case,
#'   records with missing analysis values *add* records from `dataset_ref` after the
#'   data are sorted by the variables in `by_vars` and by visit (`AVISITN` and `AVISIT`),
#'   as specified in the `order` argument.
#' - Variables other than `analysis_var` and `by_vars` that require LOCF (Last-Observation-
#'   Carried-Forward handling (in this case, `PARAMN`) are specified in the `keep_vars`
#'   argument.
#'
#' @code
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#'
#' advs <- tribble(
#'   ~STUDYID,  ~USUBJID,      ~VSSEQ, ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT,
#'   "CDISC01", "01-701-1015",      1, "PULSE",        1,    65,        0, "BASELINE",
#'   "CDISC01", "01-701-1015",      2, "DIABP",        2,    79,        0, "BASELINE",
#'   "CDISC01", "01-701-1015",      3, "DIABP",        2,    80,        2, "WEEK 2",
#'   "CDISC01", "01-701-1015",      4, "DIABP",        2,    NA,        4, "WEEK 4",
#'   "CDISC01", "01-701-1015",      5, "DIABP",        2,    NA,        6, "WEEK 6",
#'   "CDISC01", "01-701-1015",      6, "SYSBP",        3,   130,        0, "BASELINE",
#'   "CDISC01", "01-701-1015",      7, "SYSBP",        3,   132,        2, "WEEK 2"
#' )
#'
#' paramcd <- tribble(
#' ~PARAMCD,
#' "PULSE",
#' "DIABP",
#' "SYSBP"
#' )
#'
#' avisit <- tribble(
#' ~AVISITN, ~AVISIT,
#' 0, "BASELINE",
#' 2, "WEEK 2",
#' 4, "WEEK 4",
#' 6, "WEEK 6"
#' )
#'
#' advs_expected_obsv <- paramcd %>%
#' crossing(avisit)
#'
#' derive_locf_records(
#'   dataset = advs,
#'   dataset_ref = advs_expected_obsv,
#'   by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'   imputation = "add",
#'   order = exprs(AVISITN, AVISIT),
#'   keep_vars = exprs(PARAMN)
#' ) |>
#'   arrange(USUBJID, PARAMCD, AVISIT)
#'
#'
#' @caption Update records for missing analysis variable
#' @info When the `imputation` mode is set to *update*, missing `analysis_var` values
#'       are updated using values from the last record after the dataset is sorted by
#'       `by_vars` and `order`. Imputed records are added for missing timepoints (from
#'       `dataset_ref`).
#'
#' @code
#' derive_locf_records(
#'   dataset = advs,
#'   dataset_ref = advs_expected_obsv,
#'   by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'   imputation = "update",
#'   order = exprs(AVISITN, AVISIT),
#' ) |>
#'   arrange(USUBJID, PARAMCD, AVISIT)
#'
#'
#' @caption Update records for missing analysis variable while keeping the original records
#' @info When the `imputation` mode is set to *update_add*, the missing `analysis_var`
#'       values are updated using values from the last record after the dataset is sorted
#'       by `by_vars` and `order`. The updated values are added as new records, while the
#'       original records with missing `analysis_var` are retained. Imputed records are added
#'       for missing timepoints (from `dataset_ref`).
#'
#' @code
#' derive_locf_records(
#'   dataset = advs,
#'   dataset_ref = advs_expected_obsv,
#'   by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'   imputation = "update_add",
#'   order = exprs(AVISITN, AVISIT),
#' ) |>
#'   arrange(USUBJID, PARAMCD, AVISIT)
derive_locf_records <- function(dataset,
                                dataset_ref,
                                by_vars,
                                id_vars_ref = NULL,
                                analysis_var = AVAL,
                                imputation = "add",
                                order,
                                keep_vars = NULL) {
  # Input Checking
  analysis_var <- assert_symbol(enexpr(analysis_var))

  # Check if input parameters is a valid list of variables
  assert_vars(by_vars, optional = TRUE)
  assert_vars(keep_vars, optional = TRUE)
  assert_vars(id_vars_ref, optional = TRUE)
  assert_expr_list(order)

  imputation <-
    assert_character_scalar(
      imputation,
      values = c("add", "update", "update_add"),
      case_sensitive = FALSE
    )

  # Check by_vars and order variables in input datasets
  assert_data_frame(dataset_ref)
  assert_data_frame(
    dataset,
    required_vars = expr_c(
      by_vars, analysis_var, extract_vars(order), keep_vars,
      chr2vars(colnames(dataset_ref))
    )
  )


  # Setting id_vars_ref to all the variables of dataset_ref when not specified by user #
  if (is.null(id_vars_ref)) {
    id_vars_ref <- lapply(names(dataset_ref), sym)
  }


  # Prepping 'dataset_ref'
  # Get the IDs from input dataset for which the expected observations are to be added
  ids <- dataset %>%
    select(!!!setdiff(by_vars, chr2vars(colnames(dataset_ref)))) %>%
    distinct()

  exp_obsv <- ids %>%
    crossing(dataset_ref)


  # Add LOCF records
  # Get the variable names to join by
  exp_obs_by_vars <- as.character(union(by_vars, id_vars_ref))

  tmp_missing_avar <- get_new_tmp_var(exp_obsv, prefix = "tmp_missing_avar")

  # Flag the original missing analysis_var records
  dataset <- dataset %>%
    mutate(!!tmp_missing_avar := if_else(is.na(!!analysis_var), "missing", NA_character_))


  # Get unique combination of visits/timepoints per parameter per subject
  # from the input dataset
  advs_unique_original <- dataset %>%
    filter(!(is.na(!!analysis_var))) %>%
    select(all_of(exp_obs_by_vars)) %>%
    distinct()

  tmp_dtype <- get_new_tmp_var(exp_obsv, prefix = "tmp_dtype")
  tmp_new_records <- get_new_tmp_var(exp_obsv, prefix = "tmp_new_records")

  # Get the missing analysis_var (e.g., AVAL) records
  if (imputation %in% c("add", "update_add")) {
    aval_missing <- dataset %>%
      filter(is.na(!!analysis_var)) %>%
      remove_tmp_vars()
  } else {
    aval_missing <- NULL
  }

  if (imputation %in% c("update", "update_add")) {
    data_fill <- dataset

    # Get all the expected observations that are to be added to the input dataset
    exp_obsv_to_add <- exp_obsv %>%
      anti_join(dataset, by = exp_obs_by_vars) %>%
      mutate(!!tmp_new_records := "new")
  } else {
    # Get the records with missing analysis_var (e.g., AVAL) to impute
    data_fill <- dataset %>%
      filter(!is.na(!!analysis_var))

    # Get all the expected observations that are to be added to the input dataset
    exp_obsv_to_add <- exp_obsv %>%
      anti_join(advs_unique_original, by = exp_obs_by_vars) %>%
      mutate(!!tmp_new_records := "new")
  }

  # Add the expected observations to the input dataset
  # Arrange the dataset by 'order' and group it by 'by_vars'
  # Use fill() to fill the 'analysis_var' from the previous observation for the newly
  # added records

  aval_locf <- bind_rows(data_fill, exp_obsv_to_add) %>%
    mutate(!!tmp_dtype := if_else(is.na(!!analysis_var), "LOCF", NA_character_))

  if ("DTYPE" %in% colnames(aval_locf)) {
    aval_locf <- aval_locf %>%
      mutate(DTYPE = if_else(!!tmp_dtype == "LOCF", "LOCF", DTYPE, missing = DTYPE)) %>%
      select(-!!tmp_dtype)
  } else {
    aval_locf <- rename(aval_locf, DTYPE = !!tmp_dtype)
  }

  aval_locf <- aval_locf %>%
    arrange(!!!by_vars, !!!order) %>%
    group_by(!!!by_vars) %>%
    fill(!!analysis_var, !!!keep_vars) %>%
    ungroup() %>%
    filter(!(!is.na(!!tmp_missing_avar) & is.na(!!tmp_new_records) & is.na(DTYPE))) %>%
    remove_tmp_vars()


  # When imputation = 'add',  keep all variables other than analysis_var, by_vars,
  # order, id_vars_ref, keep_vars and 'DTYPE' missing. Else, keep all variables populated

  if (imputation == "add") {
    # Non-imputed records
    non_locf <- aval_locf %>%
      filter(!(DTYPE %in% c("LOCF")))

    # imputed records
    locf <- aval_locf %>%
      filter(DTYPE %in% c("LOCF"))

    aval_locf <- locf %>%
      mutate(across(
        .cols = -as.character(c(analysis_var, by_vars, order, id_vars_ref, keep_vars, "DTYPE")),
        .fns = ~ vector(typeof(.), 1)[NA]
      )) %>%
      bind_rows(non_locf)
  }

  # Output dataset:
  # If imputation == 'add' or 'update_add', add the missing analysis_var records
  # with non-missing + newly added LOCF records
  # If imputation == 'update', keep non-missing + newly added LOCF records

  bind_rows(aval_locf, aval_missing)
}
