#' Summarize Variables from an Additional Dataset Based on Conditions from Both
#' Datasets
#'
#' The function summarizes variables from an additional dataset and adds the
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
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The two datasets are joined by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
#'
#' @permitted [var_list_rename]
#'
#' @param order Sort order
#'
#'   The specified variables are used to determine the order of the records if
#'   `first_cond_lower` or `first_cond_upper` is specified or if `join_type`
#'   equals `"before"` or `"after"`.
#'
#'   If an expression is named, e.g., `exprs(EXSTDT =
#'   convert_dtc_to_dt(EXSTDTC), EXSEQ)`, a corresponding variable (`EXSTDT`) is
#'   added to the additional dataset and can be used in the filter conditions
#'   (`filter_add`, `filter_join`) and for `join_vars` and `new_vars`. The
#'   variable is not included in the output dataset.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [order_optional]
#'
#' @param new_vars Variables to add
#'
#'   The new variables can be defined by named expressions, i.e., `new_vars =
#'   exprs(<new variable> = <value>)`. The value must be defined such that it
#'   results in a single record per by group, e.g., by using a summary function
#'   like `mean()`, `sum()`, ...
#'
#' @permitted [expr_list_summary]
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   restricted additional dataset (`dataset_add` after applying `filter_add`).
#'   It is set to the observation number with respect to `order`. For each by
#'   group (`by_vars`) the observation number starts with `1`. The variable can
#'   be used in the conditions (`filter_join`, `first_cond_upper`,
#'   `first_cond_lower`). It can also be used to select consecutive observations
#'   or the last observation.
#'
#'   The variable is not included in the output dataset. To include it specify
#'   it for `new_vars`.
#'
#' @permitted [var]
#'
#' @param join_vars Variables to use from additional dataset
#'
#'   Any extra variables required from the additional dataset for `filter_join`
#'   should be specified for this argument. Variables specified for `new_vars`
#'   do not need to be repeated for `join_vars`. If a specified variable exists
#'   in both the input dataset and the additional dataset, the suffix ".join" is
#'   added to the variable from the additional dataset.
#'
#'   If an expression is named, e.g., `exprs(EXSTDT =
#'   convert_dtc_to_dt(EXSTDTC))`, a corresponding variable is added to the
#'   additional dataset and can be used in the filter conditions (`filter_add`,
#'   `filter_join`) and for `new_vars`.
#'
#'   The variables are not included in the output dataset.
#'
#' @permitted [var_expr_list]
#'
#' @param join_type Observations to keep after joining
#'
#'   The argument determines which of the joined observations are kept with
#'   respect to the original observation. For example, if `join_type = "after"`
#'   is specified all observations after the original observations are kept.
#'
#' @permitted [join_type]
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
#' @param check_type Check uniqueness?
#'
#'   If `"message"`, `"warning"` or `"error"` is specified, the specified
#'   message is issued if the observations of the input dataset (`dataset`) or
#'   the restricted additional dataset (`dataset_add` after applying
#'   `filter_add`) are not unique with respect to the by variables and the
#'   order.
#'
#'   The uniqueness is checked only if `tmp_obs_nr_var`, `first_cond_lower`,
#'   or `first_cond_upper` is specified or `join_type` equals `"before"` or
#'   `"after"`.
#'
#' @permitted [msg_type]
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
#' 1. The variables specified for `new_vars` are created and merged to the input
#' dataset. I.e., the output dataset contains all observations from the input
#' dataset. For observations without a matching observation in the joined
#' dataset the new variables are set as specified by `missing_values` (or to
#' `NA` for variables not in `missing_values`). Observations in the additional
#' dataset which have no matching observation in the input dataset are ignored.
#'
#' `r roxygen_save_memory()`
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars`
#'   derived from the additional dataset (`dataset_add`).
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
#' @info The examples focus on the functionality specific to this function. For
#'   examples of functionality common to all "joined" functions like
#'   `filter_join`, `filter_add`, `join_vars`, ... please see the examples
#'   of [derive_vars_joined()].
#'
#' @caption Derive cumulative dose before event (`CUMDOSA`)
#'
#' @info Deriving the cumulative actual dose up to the day of the adverse event
#'   in the `ADAE` dataset.
#'
#'   - `USUBJID` is specified for `by_vars` to join the `ADAE` and the `ADEX`
#'    dataset by subject.
#'   - `filter_join` is specified to restrict the `ADEX` dataset to the days up
#'    to the adverse event. `ADY.join` refers to the study day in `ADEX`.
#'   - The new variable `CUMDOSA` is defined by the `new_vars` argument. It is
#'   set to the sum of `AVAL`.
#'   - As `ADY` from `ADEX` is used in `filter_join` (but not in `new_vars`), it
#'    needs to be specified for `join_vars`.
#'    - The `join_type` is set to `"all"` to consider all records in the joined
#'    dataset. `join_type = "before"` can't by used here because then doses at
#'    the same day as the adverse event would be excluded.
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
#'   new_vars = exprs(CUMDOSA = sum(AVAL, na.rm = TRUE))
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
#' @caption Selecting records (`join_type = "before"`, `join_type = "after"`)
#'
#' @info The `join_type` argument can be used to select records from the
#'   additional dataset. For example, if `join_type = "before"` is specified,
#'   only records before the current observation are selected. If `join_type =
#'   "after"` is specified, only records after the current observation are
#'   selected.
#'
#'   To illustrate this, a variable (`SELECTED_DAYS`) is derived which contains
#'   the selected days.
#'
#' @code
#' mydata <- tribble(
#'   ~DAY,
#'   1,
#'   2,
#'   3,
#'   4,
#'   5
#' )
#'
#' derive_vars_joined_summary(
#'   mydata,
#'   dataset_add = mydata,
#'   order = exprs(DAY),
#'   join_type = "before",
#'   new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
#' )
#'
#' derive_vars_joined_summary(
#'   mydata,
#'   dataset_add = mydata,
#'   order = exprs(DAY),
#'   join_type = "after",
#'   new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
#' )
#'
#' @caption Selecting records (`first_cond_lower`, `first_cond_upper`)
#'
#' @info The `first_cond_lower` and `first_cond_upper` arguments can be used to
#'  restrict the joined dataset to a certain range of records. For example, if
#'  `first_cond_lower` is specified, the joined dataset is restricted to the
#'  last observation before the current record where the condition is
#'  fulfilled.
#'
#'  Please note:
#'  - If the condition is not fulfilled for any of the records, no records are
#'  selected.
#'  - The restriction implied by `join_type` is applied first.
#'  -  If a variable is contained in both `dataset` and `dataset_add` like `DAY`
#'  in the example below, `DAY` refers to the value from `dataset` and
#'  `DAY.join` to the value from `dataset_add`.
#'
#'  To illustrate this, a variable (`SELECTED_DAYS`) is derived which contains
#'  the selected days.
#'
#' @code
#' derive_vars_joined_summary(
#'   mydata,
#'   dataset_add = mydata,
#'   order = exprs(DAY),
#'   join_type = "before",
#'   first_cond_lower = DAY.join == 2,
#'   new_vars = exprs(SELECTED_DAYS = paste(sort(DAY), collapse = ", "))
#' )
#'
#' derive_vars_joined_summary(
#'   mydata,
#'   dataset_add = mydata,
#'   order = exprs(DAY),
#'   join_type = "after",
#'   first_cond_upper = DAY.join == 4,
#'   new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
#' )
#'
#' derive_vars_joined_summary(
#'   mydata,
#'   dataset_add = mydata,
#'   order = exprs(DAY),
#'   join_type = "all",
#'   first_cond_lower = DAY.join == 2,
#'   first_cond_upper = DAY.join == 4,
#'   new_vars = exprs(SELECTED_DAYS = paste(sort(DAY), collapse = ", "))
#' )
#'
#' @caption Derive weekly score if enough assessments are available
#'
#' @info For each planned visit the average score within the week before the
#'   visit should be derived if at least three assessments are available.
#'
#'   Please note that the condition for the number of assessments is specified
#'   in `new_vars` and not in `filter_join`. This is because the number of
#'   assessments within the week before the visit should be counted but not the
#'   number of assessments available for the subject.
#'
#' @code
#' planned_visits <- tribble(
#'   ~AVISIT,  ~ADY,
#'   "WEEK 1",    8,
#'   "WEEK 4",   29,
#'   "WEEK 8",   57
#'   ) %>%
#'   mutate(USUBJID = "1", .before = AVISIT)
#'
#' adqs <- tribble(
#'   ~ADY, ~AVAL,
#'      1,    10,
#'      2,    12,
#'      4,     9,
#'      5,     9,
#'      7,    10,
#'     25,    11,
#'     27,    10,
#'     29,    10,
#'     41,     8,
#'     42,     9,
#'     44,     5
#' ) %>%
#' mutate(USUBJID = "1")
#'
#' derive_vars_joined_summary(
#'   planned_visits,
#'   dataset_add = adqs,
#'   by_vars = exprs(USUBJID),
#'   filter_join = ADY - 7 <= ADY.join & ADY.join < ADY,
#'   join_type = "all",
#'   join_vars = exprs(ADY),
#'   new_vars = exprs(AVAL = if_else(n() >= 3, mean(AVAL, na.rm = TRUE), NA))
#' )
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

  preexisting_vars <- setdiff(chr2vars(colnames(dataset)), by_vars_left)
  duplicates <- intersect(replace_values_by_names(new_vars), preexisting_vars)
  if (length(duplicates) > 0) {
    cli_abort(c(
      paste(
        "The variable{?s} {.var {duplicates}} in {.arg new_vars} {?is/are} already",
        "in {.arg dataset}"
      ),
      "Please make appropriate modifications to {.arg dataset} or {.arg new_vars}."
    ))
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
    derive_vars_merged_summary(
      dataset_add = data_joined,
      by_vars = exprs(!!!by_vars_left, !!tmp_obs_nr),
      new_vars = add_suffix_to_vars(new_vars, vars = common_vars, suffix = ".join"),
      missing_values = missing_values
    ) %>%
    remove_tmp_vars()
}
