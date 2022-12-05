#' Filter the First or Last Observation for Each By Group
#'
#' Filters the first or last observation for each by group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
#'   expected.
#'
#' @param by_vars Grouping variables
#'
#'   *Default*: `NULL`
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
#' @author Stefan Bundfuss
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
#' library(admiral.test)
#' data("admiral_ex")
#'
#' # Select first dose for each patient
#' admiral_ex %>%
#'   filter_extreme(
#'     by_vars = vars(USUBJID),
#'     order = vars(EXSEQ),
#'     mode = "first"
#'   ) %>%
#'   select(USUBJID, EXSEQ)
#'
#' # Select highest dose for each patient on the active drug
#' admiral_ex %>%
#'   filter(EXTRT != "PLACEBO") %>%
#'   filter_extreme(
#'     by_vars = vars(USUBJID),
#'     order = vars(EXDOSE),
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
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )

  # group and sort input dataset
  tmp_obs_nr <- get_new_tmp_var(dataset)
  if (!is.null(by_vars)) {
    assert_has_variables(dataset, vars2chr(by_vars))

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
