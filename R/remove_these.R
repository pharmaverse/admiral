#' Turn a List of Quosures into a Character Vector
#'
#' @param quosures A `list` of `quosures` created using [`vars()`]
#'
#' @return A character vector
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' vars2chr(vars(USUBJID, AVAL))
vars2chr <- function(quosures) {
  rlang::set_names(
    map_chr(quosures, ~ as_string(quo_get_expr(.x))),
    names(quosures)
  )
}

#' Negate List of Variables
#'
#' The function adds a minus sign as prefix to each variable.
#'
#' This is useful if a list of variables should be removed from a dataset,
#' e.g., `select(!!!negate_vars(by_vars))` removes all by variables.
#'
#' @param vars List of variables created by `vars()`
#'
#' @return A list of `quosures`
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' negate_vars(vars(USUBJID, STUDYID))
negate_vars <- function(vars = NULL) {
  assert_vars(vars, optional = TRUE)
  if (is.null(vars)) {
    NULL
  } else {
    lapply(vars, function(var) expr(-!!quo_get_expr(var)))
  }
}

#' Optional Filter
#'
#' Filters the input dataset if the provided expression is not `NULL`
#'
#' @param dataset Input dataset
#' @param filter A filter condition. Must be a quosure.
#'
#' @return A `data.frame` containing all rows in `dataset` matching `filter` or
#' just `dataset` if `filter` is `NULL`
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @keywords user_utility
#'
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset)
  assert_filter_cond(filter, optional = TRUE)

  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}
