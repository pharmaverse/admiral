#' Add the First or Last Observation for Each By Group as New Records
#'
#' Add the first or last observation for each by group as new observations. It
#' can be used for example for adding the maximum or minimum value as a separate
#' visit. All variables of the selected observation are kept. This distinguish
#' `derive_extreme_records()` from `derive_summary_records()`, where only the by
#' variables are populated for the new records.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is added
#'   to the input dataset. If `"last"` is specified, the last observation of
#'   each by group is added to the input dataset.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param filter Filter for observations to consider
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for selecting the first or last observation. If the parameter is not
#'   specified, all observations are considered.
#'
#'   *Default*: `NULL`
#'
#'   *Permitted Values*: a condition
#'
#' @inheritParams filter_extreme
#' @inheritParams derive_summary_records
#'
#' @details
#'   1. The input dataset is restricted as specified by the `filter` parameter.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is selected.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the selected observations.
#'   1. The observations are added to input dataset.
#'
#'
#' @return The input dataset with the first or last observation of each by group
#'   added as new observations.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adlb <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
#'   "1",      1,          113,      1,
#'   "1",      2,          113,      2,
#'   "1",      3,          117,      3,
#'   "2",      1,          101,      1,
#'   "2",      2,          101,      2,
#'   "2",      3,           95,      3
#' )
#'
#' # Add a new record for each USUBJID storing the minimum value (first AVAL).
#' # If multiple records meet the minimum criterion, take the first value by
#' # AVISITN. Set AVISITN = 97 and DTYPE = MINIMUM for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL, AVISITN),
#'   mode = "first",
#'   filter = !is.na(AVAL),
#'   set_values_to = exprs(
#'     AVISITN = 97,
#'     DTYPE = "MINIMUM"
#'   )
#' )
#'
#' # Add a new record for each USUBJID storing the maximum value (last AVAL).
#' # If multiple records meet the maximum criterion, take the first value by
#' # AVISITN. Set AVISITN = 98 and DTYPE = MAXIMUM for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(desc(AVAL), AVISITN),
#'   mode = "first",
#'   filter = !is.na(AVAL),
#'   set_values_to = exprs(
#'     AVISITN = 98,
#'     DTYPE = "MAXIMUM"
#'   )
#' )
#'
#' # Add a new record for each USUBJID storing for the last value.
#' # Set AVISITN = 99 and DTYPE = LOV for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISITN = 99,
#'     DTYPE = "LOV"
#'   )
#' )
derive_extreme_records <- function(dataset = NULL,
                                   dataset_add = NULL,
                                   dataset_ref = NULL,
                                   by_vars = NULL,
                                   order = NULL,
                                   mode = NULL,
                                   filter_add = NULL,
                                   check_type = "warning",
                                   exist_flag = NULL,
                                   true_value = "Y",
                                   false_value = "N",
                                   set_values_to,
                                   filter) {
  if (!missing(filter)) {
    deprecate_warn("0.11.0", "derive_extreme_records(filter = )", "derive_extreme_records(filter_add = )")
    filter_add <- enexpr(filter)
  }

  # Check input parameters
  assert_vars(by_vars, optional = TRUE)
  assert_order_vars(order, optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = expr_c(
      by_vars, extract_vars(order)
    ),
    optional = TRUE
  )
  mode <- assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE,
    optional = TRUE
  )
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  exist_flag <- assert_symbol(enexpr(exist_flag), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_varval_list(set_values_to)
  if (is.null(dataset) && is.null(dataset_add)) {
    abort(paste(
      "Neither `dataset` nor `dataset_add` is specified.",
      "At least one of them must be specified.",
      sep = "\n")
    )
  }

  # Create new observations
  if (is.null(dataset_add)) {
    dataset_add <- dataset
  }
  new_add_obs <- filter_if(dataset_add, filter_add)

  if (!is.null(order)) {
    new_add_obs <- filter_extreme(
      new_add_obs,
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    )
  }

  if (!is.null(dataset_ref)) {
    add_vars <- colnames(dataset_add)
    ref_vars <- colnames(dataset_ref)

    new_ref_obs <- anti_join(
      select(dataset_ref, intersect(add_vars, ref_vars)),
      select(new_add_obs, !!!by_vars),
      by = sapply(by_vars, as_name) # nolint: undesirable_function_linter
    )

    if (!is.null(exist_flag)) {
      new_add_obs <- mutate(new_add_obs, !!exist_flag := true_value)
      new_ref_obs <- mutate(new_ref_obs, !!exist_flag := false_value)
    }
    new_obs <- bind_rows(new_add_obs, new_ref_obs)
  } else {
    new_obs <- new_add_obs
  }

  new_obs <- process_set_values_to(
    new_obs,
    set_values_to = set_values_to
  )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
