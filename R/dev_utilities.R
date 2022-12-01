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
#' @author Thomas Neitmann
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
#' @author Ondrej Slama
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
#' @author Thomas Neitmann, Ondrej Slama
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
    (expr[[1L]] == quote(enquo) || expr[[1L]] == quote(rlang::enquo)) &&
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

#' Extract All Symbols from a List of Quosures
#'
#' @param x An `R` object
#' @param side One of `"lhs"` (the default) or `"rhs"`
#'
#' @return A list of `quosures`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#' @family dev_utility
#' @export
extract_vars <- function(x, side = "lhs") {
  if (is.null(x)) {
    NULL
  } else if (is.list(x)) {
    do.call(quo_c, map(x, extract_vars, side))
  } else if (is_quosure(x)) {
    env <- quo_get_env(x)
    symbols <- syms(all.vars(quo_get_expr(x)))
    map(symbols, ~ quo_set_env(quo(!!.x), env))
  } else if (is_formula(x)) {
    funs <- list("lhs" = f_lhs, "rhs" = f_rhs)
    assert_character_scalar(side, values = names(funs))
    quo_set_env(
      quo(!!funs[[side]](x)),
      env = attr(x, ".Environment")
    )
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


#' Turn a Quosure into a String
#'
#' @details
#' This function is missing in earlier version of {rlang} which is why we re-
#' implement it here.
#'
#' @param x A `quosure`
#'
#' @return A `character` vector
#'
#' @keywords dev_utility
#' @family dev_utility
#'
#' @export
as_name <- function(x) {
  if (is_quosure(x)) {
    x <- quo_get_expr(x)
  }
  as_string(x)
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

#' check that argument contains valid variable(s) created with `vars()` or
#' Source Variables from a List of Quosures
#'
#' @param arg A function argument to be checked
#'
#' @return A TRUE if variables were valid variable
#'
#' @export
#'
#' @keywords dev_utility
#' @family dev_utility
contains_vars <- function(arg) {
  inherits(arg, "quosures") && all(map_lgl(arg, quo_is_symbol) | names(arg) != "")
}

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
#' @keywords dev_utility
#' @family dev_utility
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' vars2chr(vars(USUBJID, AVAL))
vars2chr <- function(quosures) {
  rlang::set_names(
    map_chr(quosures, ~ as_string(quo_get_expr(.x))),
    names(quosures)
  )
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
#' @keywords dev_utility
#' @family dev_utility
#'
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_filter_cond(filter, optional = TRUE)
  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}
