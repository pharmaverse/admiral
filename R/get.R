#' Get Constant Variables
#'
#' @param dataset A data frame.
#' @param by_vars By variables
#'   The groups defined by the by variables are considered separately. I.e., if
#'   a variable is constant within each by group, it is returned.
#'
#' @param ignore_vars Variables to ignore
#'   The specified variables are not considered, i.e., they are not returned
#'   even if they are constant (unless they are included in the by variables).
#'
#'   *Permitted Values:* A list of variable names or selector function calls
#'   like `starts_with("EX")`
#'
#' @keywords get
#' @family get
#'
#'
#' @return Variable vector.
#' @export
get_constant_vars <- function(dataset, by_vars, ignore_vars = NULL) {
  assert_data_frame(dataset, optional = FALSE)
  assert_vars(by_vars, optional = FALSE)
  assert_vars(ignore_vars, optional = TRUE)

  non_by_vars <- setdiff(names(dataset), vars2chr(by_vars))

  if (!is.null(ignore_vars)) {
    non_by_vars <- setdiff(
      non_by_vars,
      vars_select(non_by_vars, !!!ignore_vars)
    )
  }

  # get unique values within each group by variables
  unique_count <- dataset %>%
    group_by(!!!by_vars) %>%
    summarise(across(!!non_by_vars, n_distinct)) %>%
    ungroup() %>%
    select(!!!syms(non_by_vars))

  # determine variables which are constant within each by group
  constant_vars <- unique_count %>%
    map_lgl(~ all(.x == 1)) %>%
    which() %>%
    names() %>%
    syms()

  exprs(!!!by_vars, !!!constant_vars)
}


#' Get Duplicates From a Vector
#'
#' @param x An atomic vector
#'
#' @return A vector of the same type as `x` contain duplicate values
#' @export
#'
#' @keywords get
#' @family get
#'
#' @export
#'
#' @examples
#' get_duplicates(1:10)
#'
#' get_duplicates(c("a", "a", "b", "c", "d", "d"))
get_duplicates <- function(x) {
  assert_atomic_vector(x)

  unique(x[duplicated(x)])
}

#' Get Source Variables from a List of Expressions
#'
#' @param expressions A list of expressions
#'
#' @param quosures *Deprecated*, please use `expressions` instead.
#'
#'
#' @keywords get
#' @family get
#'
#' @return A list of expressions
#' @export
get_source_vars <- function(expressions, quosures) {
  if (!missing(quosures)) {
    deprecate_warn(
      "0.10.0",
      "get_source_vars(quosures = )",
      "get_source_vars(expressions = )"
    )
    expressions <- map(quosures, rlang::quo_get_expr)
  }
  assert_varval_list(expressions, optional = TRUE)

  source_vars <- expr_c(expressions)[lapply(expr_c(expressions), is.symbol) == TRUE]

  if (length(source_vars) == 0) {
    NULL
  } else {
    source_vars
  }
}
