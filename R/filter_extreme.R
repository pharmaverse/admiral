#' Filter out the first or last observation for each group
#'
#' Filters out the first or last observation for each group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
#'   expected.
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is
#'   included in the output dataset. If `"last"` is specified, the last
#'   observation of each by group is included in the output dataset.
#'
#'   Permitted Values:  `"first"`, `"last"`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   Default: `"none"`
#'
#'   Permitted Values: `"none"`, `"warning"`, `"error"`
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
#' @keywords adam
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(magrittr)
#'
#' data("ex")
#'
#' # selecting first dose for each patient
#' filter_extreme(ex,
#'                order = exprs(EXSEQ),
#'                by_vars = exprs(USUBJID),
#'                mode = 'first') %>%
#'   select(USUBJID, EXSEQ)
#'
#' # selecting highest dose for each patient
#' filter_extreme(ex,
#'                order = exprs(EXDOSE),
#'                by_vars = exprs(USUBJID),
#'                check_type = "none") %>%
#'   select(USUBJID, EXDOSE)
#'

filter_extreme <- function(dataset,
                           order,
                           by_vars,
                           mode = "last",
                           check_type = "warning") {
  arg_match(mode, c("first", "last"))
  arg_match(check_type, c("none", "warning", "error"))

  # group and sort input dataset
  if (!missing(by_vars)) {
    assert_has_variables(dataset, map_chr(by_vars, as_string))

    data <- dataset %>%
      derive_obs_number(order = order,
                        by_vars = by_vars,
                        check_type = check_type) %>%
      group_by(!!!by_vars)
  }
  else {
    data <- dataset %>% derive_obs_number(order = order,
                                          check_type = check_type)
  }

  if (mode == "first") {
    # select first observation (for each group)
    data <- data %>%
      slice(1)
  }
  else {
    # select last observation (for each group)
    data <- data %>%
      slice(n())
  }
  data %>%
    ungroup() %>%
    select(-starts_with("temp_"))
}
