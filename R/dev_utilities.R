#' Negated Value Matching
#'
#' Returns a `logical` vector indicating if there is *no* match of the
#' left operand in the right operand.
#'
#' @param x The values to be matched
#' @param table The values to be matched against
#'
#' @return A `logical` vector
#'
#'
#' @keywords dev_utility
#' @family dev_utility
#' @export
`%notin%` <- function(x, table) { # nolint
  !(x %in% table)
}

#' Helper Function to Convert Date (or Date-time) Objects to Characters of dtc Format
#' (-DTC type of variable)
#'
#' @param dtm date or date-time
#'
#' @return `character` vector
#'
#'
#' @keywords dev_utility
#' @family dev_utility
#' @export
convert_dtm_to_dtc <- function(dtm) {
  stopifnot(lubridate::is.instant(dtm))
  format(dtm, "%Y-%m-%dT%H:%M:%S")
}

#' Extract Argument Name from an Expression
#'
#' @param expr An expression created inside a function using `substitute()`
#'
#'
#' @return `character` vector
#'
#' @keywords dev_utility
#' @family dev_utility
#'
#' @export
arg_name <- function(expr) { # nolint
  if (length(expr) == 1L && is.symbol(expr)) {
    deparse(expr)
  } else if (length(expr) == 2L &&
    (expr[[1L]] == quote(enexpr) || expr[[1L]] == quote(rlang::enexpr)) &&
    is.symbol(expr[[2L]])) {
    deparse(expr[[2L]])
  } else if (is.call(expr) && length(expr) >= 2 && is.symbol(expr[[2]])) {
    deparse(expr[[2L]])
  } else if (is.call(expr) && length(expr) >= 2 && is.call(expr[[2]])) {
    arg_name(expr[[2L]])
  } else {
    abort(paste0("Could not extract argument name from `", deparse(expr), "`"))
  }
}

#' Extract All Symbols from a List of Expressions
#'
#' @param x An `R` object
#' @param side One of `"lhs"` (the default) or `"rhs"`
#'
#' @return A list of expressions
#'
#'
#' @keywords dev_utility
#' @family dev_utility
#' @export
extract_vars <- function(x, side = "lhs") {
  if (is.null(x)) {
    NULL
  } else if (is.list(x)) {
    do.call(expr_c, map(x, extract_vars, side))
  } else if (is_expression(x)) {
    syms(all.vars(x))
  } else if (is_formula(x)) {
    funs <- list("lhs" = f_lhs, "rhs" = f_rhs)
    assert_character_scalar(side, values = names(funs))
    expr(!!funs[[side]](x))
  } else {
    abort()
  }
}

#' Or
#'
#' @param lhs Any valid R expression
#' @param rhs Any valid R expression
#'
#' @details
#' The function evaluates the expression `lhs` and if this expression results
#' in an error, it catches that error and proceeds with evaluating the expression
#' `rhs` and returns that result.
#'
#' @return Either the result of evaluating `lhs`, `rhs` or an error
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
`%or%` <- function(lhs, rhs) {
  tryCatch(lhs, error = function(e) rhs)
}

#' Valid Time Units
#'
#' Contains the acceptable character vector of valid time units
#'
#' @return A `character` vector of valid time units
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
valid_time_units <- function() {
  c("years", "months", "days", "hours", "minutes", "seconds")
}

#' check that argument contains valid variable(s) created with `exprs()` or
#' Source Variables from a List of Expressions
#'
#' @param arg A function argument to be checked
#'
#' @return A `TRUE` if variables were valid variable
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
contains_vars <- function(arg) {
  inherits(arg, "list") && all(map_lgl(arg, is_symbol) | names(arg) != "")
}

#' Turn a List of Expressions into a Character Vector
#'
#' @param expressions A `list` of expressions created using [`exprs()`]
#'
#' @param quosures *Deprecated*, please use `expressions` instead.
#'
#' @return A character vector
#'
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#'
#' vars2chr(exprs(USUBJID, AVAL))
vars2chr <- function(expressions, quosures) {
  if (!missing(quosures)) {
    deprecate_warn(
      "0.10.0",
      "vars2chr(quosures = )",
      "vars2chr(expressions = )"
    )
    expressions <- map(quosures, rlang::quo_get_expr)
  }
  rlang::set_names(
    map_chr(expressions, as_string),
    names(expressions)
  )
}

#' Optional Filter
#'
#' Filters the input dataset if the provided expression is not `NULL`
#'
#' @param dataset Input dataset
#' @param filter A filter condition. Must be an expression.
#'
#' @return A `data.frame` containing all rows in `dataset` matching `filter` or
#' just `dataset` if `filter` is `NULL`
#'
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
#'
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_filter_cond(filter, optional = TRUE)
  if (is.null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}
