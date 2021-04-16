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
#' data("ex")
#'
#' # selecting first dose for each patient
#' filter_extreme(ex,
#'                order = rlang::exprs(EXSEQ),
#'                by_vars = rlang::exprs(USUBJID),
#'                mode = 'first')
#'
#' # selecting highest dose for each patient
#' filter_extreme(ex,
#'                order = rlang::exprs(EXDOSE),
#'                by_vars = rlang::exprs(USUBJID))
#'

filter_extreme <- function(dataset,
                           order,
                           by_vars,
                           mode = "last"){
  arg_match(mode, c("first", "last"))

  # group and sort input dataset
  if (!missing(by_vars)){
    assert_has_variables(dataset, map_chr(by_vars, as_string))

    data <- dataset %>%
      derive_obs_number(order = order,
                        by_vars = by_vars,
                        check_type = "warning") %>%
      group_by(!!!by_vars)
  }
  else{
    data <- dataset %>% derive_obs_number(order = order,
                                          check_type = "warning")
  }

  if (mode == "first"){
    # select first observation (for each group)
    data <- data %>%
      slice(1)
  }
  else{
    # select last observation (for each group)
    data <- data %>%
      slice(n())
  }
  data %>%
    ungroup() %>%
    select(-starts_with("temp_"))
}
