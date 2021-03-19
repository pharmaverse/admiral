#' Filter out the first observation for each group
#'
#' Filters out the first observation for each group.
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
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#'   Permitted Values: list of variables
#'
#' @author Stefan Bundfuss
#'
#' @return A dataset containing the first observation of each by group
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' data("ex")
#'
#' # selecting first dose for each patient
#' filter_first(ex,
#'              order = exprs(EXSEQ),
#'              by_vars = exprs(USUBJID))
#'
#' # selecting highest dose for each patient
#' filter_first(ex,
#'              order = exprs(desc(EXDOSE)),
#'              by_vars = exprs(USUBJID))
#'

filter_first <- function(dataset,
                         order,
                         by_vars){
  if (!missing(by_vars)){
    data <- dataset %>% group_by(!!!by_vars) %>%
      arrange(!!!order, .by_group = TRUE)
  }
  else{
    data <- dataset %>%
      arrange(!!!order)
  }
  data %>%
    slice(1)
}
