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
#' @author Stefan Bundfuss
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
#'   by_vars = vars(USUBJID),
#'   order = vars(AVAL, AVISITN),
#'   mode = "first",
#'   filter = !is.na(AVAL),
#'   set_values_to = vars(
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
#'   by_vars = vars(USUBJID),
#'   order = vars(desc(AVAL), AVISITN),
#'   mode = "first",
#'   filter = !is.na(AVAL),
#'   set_values_to = vars(
#'     AVISITN = 98,
#'     DTYPE = "MAXIMUM"
#'   )
#' )
#'
#' # Add a new record for each USUBJID storing for the last value.
#' # Set AVISITN = 99 and DTYPE = LOV for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = vars(USUBJID),
#'   order = vars(AVISITN),
#'   mode = "last",
#'   set_values_to = vars(
#'     AVISITN = 99,
#'     DTYPE = "LOV"
#'   )
#' )
derive_extreme_records <- function(dataset,
                                   by_vars = NULL,
                                   order,
                                   mode,
                                   check_type = "warning",
                                   filter = NULL,
                                   set_values_to) {
  # Check input parameters
  assert_vars(by_vars, optional = TRUE)
  assert_order_vars(order)
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_varval_list(set_values_to)

  # Create new observations
  new_obs <- dataset %>%
    filter_if(filter) %>%
    filter_extreme(
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    ) %>%
    mutate(!!!set_values_to)

  # Create output dataset
  bind_rows(dataset, new_obs)
}
