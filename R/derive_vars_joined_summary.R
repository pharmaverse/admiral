#' Summarize Variables from an Additional Dataset Based on Conditions from Both
#' Datasets
#'
#' The function summarises variables from an additional dataset and adds the
#' summarized values as new variables to the input dataset. The selection of the
#' observations from the additional dataset can depend on variables from both
#' datasets. For example, all doses before the current observation can be
#' selected and the sum be added to the input dataset.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, the `join_vars`,
#'   and the `order` argument are expected.
#'
#' @param by_vars Grouping variables
#'
#'   The two datasets are joined by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each observation of the
#'   input dataset the first or last observation from the joined dataset is
#'   selected with respect to the specified order. The specified variables are
#'   expected in the additional dataset (`dataset_add`). If a variable is
#'   available in both `dataset` and `dataset_add`, the one from `dataset_add`
#'   is used for the sorting.
#'
#'   If an expression is named, e.g., `exprs(EXSTDT =
#'   convert_dtc_to_dt(EXSTDTC), EXSEQ)`, a corresponding variable (`EXSTDT`) is
#'   added to the additional dataset and can be used in the filter conditions
#'   (`filter_add`, `filter_join`) and for `join_vars` and `new_vars`. The
#'   variable is not included in the output dataset.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted list of expressions created by `exprs()`, e.g.,
#'    `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @param new_vars Variables to add
#'
#'   The specified variables from the additional dataset are added to the output
#'   dataset. Variables can be renamed by naming the element, i.e., `new_vars =
#'   exprs(<new name> = <old name>)`.
#'
#'   For example `new_vars = exprs(var1, var2)` adds variables `var1` and `var2`
#'   from `dataset_add` to the input dataset.
#'
#'   And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
#'   `old_var2` from `dataset_add` and adds them to the input dataset renaming
#'   `old_var2` to `new_var2`.
#'
#'   Values of the added variables can be modified by specifying an expression.
#'   For example, `new_vars = LASTRSP = exprs(str_to_upper(AVALC))` adds the
#'   variable `LASTRSP` to the dataset and sets it to the upper case value of
#'   `AVALC`.
#'
#'   If the argument is not specified or set to `NULL`, all variables from the
#'   additional dataset (`dataset_add`) are added.
#'
#' @permitted list of variables or named expressions created by `exprs()`
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   additional dataset (`dataset_add`). It is set to the observation number
#'   with respect to `order`. For each by group (`by_vars`) the observation
#'   number starts with `1`. The variable can be used in the conditions
#'   (`filter_join`, `first_cond_upper`, `first_cond_lower`). It can also be
#'   used to select consecutive observations or the last observation.
#'
#'   The variable is not included in the output dataset. To include it specify
#'   it for `new_vars`.
#'
#' @param join_vars Variables to use from additional dataset
#'
#'   Any extra variables required from the additional dataset for `filter_join`
#'   should be specified for this argument. Variables specified for `new_vars`
#'   do not need to be repeated for `join_vars`. If a specified variable exists
#'   in both the input dataset and the additional dataset, the suffix ".join" is
#'   added to the variable from the additional dataset.
#'
#'   If an expression is named, e.g., `exprs(EXTDT =
#'   convert_dtc_to_dt(EXSTDTC))`, a corresponding variable is added to the
#'   additional dataset and can be used in the filter conditions (`filter_add`,
#'   `filter_join`) and for `new_vars`. The variable is not included in the
#'   output dataset.
#'
#'   The variables are not included in the output dataset.
#'
#' @permitted list of variables or named expressions created by `exprs()`
#'
#' @param first_cond_lower Condition for selecting range of data (before)
#'
#'   If this argument is specified, the other observations are restricted from
#'   the first observation before the current observation where the specified
#'   condition is fulfilled up to the current observation. If the condition is
#'   not fulfilled for any of the other observations, no observations are
#'   considered.
#'
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only from a
#'   certain observation before the current observation up to the current
#'   observation. For an example see the last example below.
#'
#' @permitted [condition]
#'
#' @param first_cond_upper Condition for selecting range of data (after)
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the other observations, no
#'   observations are considered.
#'
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment. For an example see the last example below.
#'
#' @permitted [condition]
#'
#' @param filter_join Filter for the joined dataset
#'
#'   The specified condition is applied to the joined dataset. Therefore
#'   variables from both datasets `dataset` and `dataset_add` can be used.
#'
#'   Variables created by `order` or `new_vars` arguments can be used in the
#'   condition.
#'
#'   The condition can include summary functions like `all()` or `any()`. The
#'   joined dataset is grouped by the original observations.
#'
#' @permitted [condition]
#'
#' @inheritParams get_joined_data
#' @inheritParams derive_vars_merged
#'
#' @details
#'
#' 1. The variables specified by `order` are added to the additional dataset
#' (`dataset_add`).
#'
#' 1. The variables specified by `join_vars` are added to the additional dataset
#' (`dataset_add`).
#'
#' 1. The records from the additional dataset (`dataset_add`) are restricted to
#' those matching the `filter_add` condition.
#'
#' 1. The input dataset and the (restricted) additional dataset are left joined
#' by the grouping variables (`by_vars`). If no grouping variables are
#' specified, a full join is performed.
#'
#' 1. If `first_cond_lower` is specified, for each observation of the input
#'     dataset the joined dataset is restricted to observations from the first
#'     observation where `first_cond_lower` is fulfilled (the observation fulfilling
#'     the condition is included) up to the observation of the input dataset. If for
#'     an observation of the input dataset the condition is not fulfilled, the
#'     observation is removed.
#'
#'     If `first_cond_upper` is specified, for each observation of the input
#'     dataset the joined dataset is restricted to observations up to the first
#'     observation where `first_cond_upper` is fulfilled (the observation
#'     fulfilling the condition is included). If for an observation of the input
#'     dataset the condition is not fulfilled, the observation is removed.
#'
#'     For an example see the last example in the "Examples" section.
#'
#' 1. The joined dataset is restricted by the `filter_join` condition.
#'
#' 1. If `order` is specified, for each observation of the input dataset the
#' first or last observation (depending on `mode`) is selected.
#'
#' 1. The variables specified for `new_vars` are created (if requested) and
#' merged to the input dataset. I.e., the output dataset contains all
#' observations from the input dataset. For observations without a matching
#' observation in the joined dataset the new variables are set as specified by
#' `missing_values` (or to `NA` for variables not in `missing_values`).
#' Observations in the additional dataset which have no matching observation in
#' the input dataset are ignored.
#'
#' `r roxygen_save_memory()`
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars` from
#'   the additional dataset (`dataset_add`).
#'
#' @seealso [derive_vars_joined()], [derive_vars_merged_summary()],
#' [derive_var_joined_exist_flag()], [filter_joined()]
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @export
#'
#' @examplesx
#'
#' @caption Derive cumulative dose before event (`CUMDOSA`)
#'
#' @info To the `ADAE` dataset the cumulative actual dose up to the day of the
#'   adverse event should be added.
#'
#'   - `USUBJID` is specified for `by_vars` to join the `ADAE` and the `ADEX`
#'    dataset by subject.
#'   - `filter_join` is specified to restrict the `ADEX` dataset to the days up
#'    to the adverse event. `ADY.join` refers to the study day in `ADEX`.
#'   - The new variable `CUMDOSE` is defined by the `new_vars` argument. It is
#'   set to the sum of `AVAL`.
#'   - As `ADY` from `ADEX` is used in `filter_join` (but not in `new_vars`), it
#'    needs to be specified for `join_vars`.
#'
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adex <- tribble(
#'   ~USUBJID, ~ADY, ~AVAL,
#'   "1",         1,    10,
#'   "1",         8,    20,
#'   "1",        15,    10,
#'   "2",         8,     5
#' )
#'
#' adae <- tribble(
#'   ~USUBJID, ~ADY, ~AEDECOD,
#'   "1",         2, "Fatigue",
#'   "1",         9, "Influenza",
#'   "1",        15, "Theft",
#'   "1",        15, "Fatigue",
#'   "2",         4, "Parasomnia",
#'   "3",         2, "Truancy"
#' )
#'
#' derive_vars_joined_summary(
#'   dataset = adae,
#'   dataset_add = adex,
#'   by_vars = exprs(USUBJID),
#'   filter_join = ADY.join <= ADY,
#'   join_type = "all",
#'   join_vars = exprs(ADY),
#'   new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE))
#' )
#'
#' @caption Define values for records without records in the additional dataset (`missing_values`)
#' @info By default, the new variables are set to `NA` for records without
#'   matching records in the restricted additional dataset. This can be changed
#'   by specifying the `missing_values` argument.
#'
#' @code
#' derive_vars_joined_summary(
#'   dataset = adae,
#'   dataset_add = adex,
#'   by_vars = exprs(USUBJID),
#'   filter_join = ADY.join <= ADY,
#'   join_type = "all",
#'   join_vars = exprs(ADY),
#'   new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE)),
#'   missing_values = exprs(CUMDOSE = 0)
#' )
#'
derive_vars_joined_summary <- function(dataset,
                               dataset_add,
                               by_vars = NULL,
                               order = NULL,
                               new_vars,
                               tmp_obs_nr_var = NULL,
                               join_vars = NULL,
                               join_type,
                               filter_add = NULL,
                               first_cond_lower = NULL,
                               first_cond_upper = NULL,
                               filter_join = NULL,
                               missing_values = NULL,
                               check_type = "warning") {
  assert_vars(by_vars, optional = TRUE)
  by_vars_left <- replace_values_by_names(by_vars)
  assert_expr_list(order, optional = TRUE)
  assert_expr_list(new_vars)
  assert_expr_list(join_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars_left)
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(
      by_vars,
      extract_vars(order),
      setdiff(extract_vars(join_vars), replace_values_by_names(order))
    )
  )

  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  first_cond_lower <- assert_filter_cond(enexpr(first_cond_lower), optional = TRUE)
  first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join), optional = TRUE)

  preexisting_vars <- chr2vars(colnames(dataset))
  preexisting_vars_no_by_vars <- preexisting_vars[which(!(preexisting_vars %in% by_vars))]
  duplicates <- intersect(replace_values_by_names(new_vars), preexisting_vars_no_by_vars)
  if (length(duplicates) > 0) {
    cli_abort(
      paste(
        "The variable{?s} {.var {duplicates}} in {.arg dataset_add} ha{?s/ve} naming",
        "conflicts with {.arg dataset}, please make the appropriate modifications",
        "to {.arg new_vars}."
      )
    )
  }

  # number observations of the input dataset to get a unique key
  # (by_vars and tmp_obs_nr)
  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr,
      by_vars = by_vars_left,
      check_type = "none"
    )

  data_joined <- get_joined_data(
    data,
    dataset_add = dataset_add,
    by_vars = by_vars,
    join_vars = expr_c(
      join_vars,
      intersect(unname(extract_vars(new_vars)), chr2vars(colnames(dataset_add)))
    ),
    join_type = join_type,
    first_cond_lower = !!first_cond_lower,
    first_cond_upper = !!first_cond_upper,
    order = order,
    tmp_obs_nr_var = !!tmp_obs_nr_var,
    filter_add = !!filter_add,
    filter_join = !!filter_join,
    check_type = check_type
  )

  common_vars <-
    chr2vars(setdiff(intersect(colnames(data), colnames(dataset_add)), vars2chr(by_vars)))
  # merge new variables to the input dataset and rename them
  data %>%
    derive_var_merged_summary(
      dataset_add = data_joined,
      by_vars = exprs(!!!by_vars_left, !!tmp_obs_nr),
      new_vars = add_suffix_to_vars(new_vars, vars = common_vars, suffix = ".join"),
      missing_values = missing_values
    ) %>%
    remove_tmp_vars()
}
