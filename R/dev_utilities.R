

#' Negated Value Matching
#'
#' Returns a `logical` vector indicating if there is *no* match of the
#' left operand in the right operand.
#'
#' @param x The values to be matched
#' @param table The values to be matched against
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_notin
#'
#' @examples
#' `%notin%` <- admiral:::`%notin%`
#' "a" %notin% c("b", "v", "k")
`%notin%` <- function(x, table) { # nolint
  !(x %in% table)
}

#' Helper Function to Convert Date (or Date-time) Objects to Characters of dtc Format
#' (-DTC type of variable)
#'
#' @param dtm date or date-time
#'
#' @return character
#'
#' @author Ondrej Slama
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_convert_dtm_to_dtc
#'
#' @examples
#' admiral:::convert_dtm_to_dtc(as.POSIXct(Sys.time()))
#' admiral:::convert_dtm_to_dtc(as.Date(Sys.time()))
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
#' @keywords dev_utility
#'
#' @rdname dev_util_arg_name
#'
#' @examples
#' test_fun <- function(something) {
#'   admiral:::arg_name(substitute(something))
#' }
#'
#' inner_function <- function(x) x
#' test_fun2 <- function(something) {
#'   admiral:::arg_name(substitute(inner_function(something)))
#' }
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
#'
#' @rdname dev_util_extract_vars
#'
#' @examples
#' admiral:::extract_vars(vars(STUDYID, USUBJID, desc(ADTM)))
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


`%or%` <- function(lhs, rhs) {
  tryCatch(lhs, error = function(e) rhs)
}

#' Replace Quosure Value with Name
#'
#' @param quosures A list of quosures
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_replace_values_by_names
#'
#' @return A list of quosures
#'
#' @examples
#' admiral:::replace_values_by_names(vars(USUBJID, TEST = VSTESTCD))
replace_values_by_names <- function(quosures) {
  vars <- map2(quosures, names(quosures), function(q, n) {
    if (n == "") {
      return(q)
    }
    quo_set_env(
      quo(!!as.symbol(n)),
      quo_get_env(q)
    )
  })
  structure(vars, class = "quosures", names = NULL)
}


#' Turn a Quosure into a String
#'
#' @details
#' This function is missing in earlier version of {rlang} which is why we re-
#' implment it here.
#'
#' @noRd
as_name <- function(x) {
  if (is_quosure(x)) {
    x <- quo_get_expr(x)
  }
  as_string(x)
}

valid_time_units <- function() {
  c("years", "months", "days", "hours", "minutes", "seconds")
}

contains_vars <- function(arg) {
  inherits(arg, "quosures") && all(map_lgl(arg, quo_is_symbol) | names(arg) != "")
}
