#' Add the First or Last Observation for Each By Group as New Records
#'
#' Adds the first or last observation for each by group as new observations. It
#' can be used for example for adding the maximum or minimum value as a separate
#' visit.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
#'   expected.
#'
#' @param by_vars Grouping variables
#'
#'   *Permitted Values:* list of variables created by `vars()`
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   *Permitted Values:* list of variables or `desc(<variable>)` function calls
#'   created by `vars()`, e.g., `vars(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is added
#'   to the input dataset. If `"last"` is specified, the last observation of
#'   each by group is added to the input dataset.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the filtered input dataset are not unique with
#'   respect to the by variables and the order.
#'
#'   *Default:* `"warning"`
#'
#'   *Permitted Values:* `"none"`, `"warning"`, `"error"`
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
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations. For example `vars(DTYPE = "MINIMUM")` defines the derivation
#'   type for the new observations.
#'
#'   *Permitted Values:* List of variable-value pairs
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
#' @keywords derivation bds
#'
#' @export
#'
#' @examples
#' adlb <- tibble::tribble(
#'   ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
#'   "1",      1,        113,   1,
#'   "1",      2,        113,   2,
#'   "1",      3,        117,   3,
#'   "2",      1,        101,   1,
#'   "2",      2,        101,   2,
#'   "2",      3,         95,   3
#' )
#'
#' # Add a new visit for the minium value (first observation if not unique)
#' derive_extreme_records(
#'   adlb,
#'   by_vars = vars(USUBJID),
#'   order = vars(AVAL, AVISITN),
#'   mode = "first",
#'   set_values_to = vars(
#'     AVISITN = 97,
#'     DTYPE = "MINIMUM"
#'   )
#' )
#'
#' # Add a new visit for the maxium value (first observation used if not unique)
#' derive_extreme_records(
#'   adlb,
#'   by_vars = vars(USUBJID),
#'   order = vars(desc(AVAL), AVISITN),
#'   mode = "first",
#'   set_values_to = vars(
#'     AVISITN = 98,
#'     DTYPE = "MAXIMUM"
#'   )
#' )
#'
#' # Add a new visit for the last value
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
    filter_extreme(by_vars = by_vars,
                   order = order,
                   mode = mode,
                   check_type = check_type) %>%
    mutate(!!!set_values_to)

  # Create output dataset
  bind_rows(dataset, new_obs)


}
