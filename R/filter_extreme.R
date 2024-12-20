#' Filter the First or Last Observation for Each By Group
#'
#' Filters the first or last observation for each by group.
#'
#' @param dataset `r roxygen_param_dataset(expected_vars = c("by_vars", "order"))`
#'
#' @param by_vars Grouping variables
#'
#'   *Default*: `NULL`
#'
#'   `r roxygen_param_by_vars()`
#'
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is
#'   included in the output dataset. If `"last"` is specified, the last
#'   observation of each by group is included in the output dataset.
#'
#'   *Permitted Values:*  `"first"`, `"last"`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   *Default:* `"warning"`
#'
#'   *Permitted Values:* `"none"`, `"warning"`, `"error"`
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is included in the output dataset.
#'
#'
#' @return A dataset containing the first or last observation of each by group
#'
#' @family utils_fil
#'
#' @keywords utils_fil
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' ex <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~EXSEQ, ~EXDOSE,    ~EXTRT,
#'   "PILOT01",    "EX", "01-1442",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "01-1442",      2,      54,    "XANO",
#'   "PILOT01",    "EX", "01-1442",      3,      54,    "XANO",
#'   "PILOT01",    "EX", "01-1444",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "01-1444",      2,      81,    "XANO",
#'   "PILOT01",    "EX", "05-1382",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "08-1213",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "10-1053",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "10-1053",      2,      54,    "XANO",
#'   "PILOT01",    "EX", "10-1183",      1,       0, "PLACEBO",
#'   "PILOT01",    "EX", "10-1183",      2,       0, "PLACEBO",
#'   "PILOT01",    "EX", "10-1183",      3,       0, "PLACEBO",
#'   "PILOT01",    "EX", "11-1036",      1,       0, "PLACEBO",
#'   "PILOT01",    "EX", "11-1036",      2,       0, "PLACEBO",
#'   "PILOT01",    "EX", "11-1036",      3,       0, "PLACEBO",
#'   "PILOT01",    "EX", "14-1425",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "15-1319",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "15-1319",      2,      81,    "XANO",
#'   "PILOT01",    "EX", "16-1151",      1,      54,    "XANO",
#'   "PILOT01",    "EX", "16-1151",      2,      54,    "XANO"
#' )
#'
#'
#' # Select first dose for each patient
#' ex %>%
#'   filter_extreme(
#'     by_vars = exprs(USUBJID),
#'     order = exprs(EXSEQ),
#'     mode = "first"
#'   ) %>%
#'   select(USUBJID, EXSEQ)
#'
#' # Select highest dose for each patient on the active drug
#' ex %>%
#'   filter(EXTRT != "PLACEBO") %>%
#'   filter_extreme(
#'     by_vars = exprs(USUBJID),
#'     order = exprs(EXDOSE),
#'     mode = "last",
#'     check_type = "none"
#'   ) %>%
#'   select(USUBJID, EXTRT, EXDOSE)
filter_extreme <- function(dataset,
                           by_vars = NULL,
                           order,
                           mode,
                           check_type = "warning") {
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error", "message"),
      case_sensitive = FALSE
    )
  assert_data_frame(dataset, required_vars = by_vars)

  # group and sort input dataset
  tmp_obs_nr <- get_new_tmp_var(dataset)
  if (!is.null(by_vars)) {
    data <- dataset %>%
      derive_var_obs_number(
        new_var = !!tmp_obs_nr,
        order = order,
        by_vars = by_vars,
        check_type = check_type
      ) %>%
      group_by(!!!by_vars)
  } else {
    data <- dataset %>%
      derive_var_obs_number(
        new_var = !!tmp_obs_nr,
        order = order,
        check_type = check_type
      )
  }

  if (mode == "first") {
    # select first observation (for each group)
    data <- data %>%
      slice(1L)
  } else {
    # select last observation (for each group)
    data <- data %>%
      slice(n())
  }
  data %>%
    ungroup() %>%
    remove_tmp_vars()
}
