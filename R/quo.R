#' Concatenate One or More Expressions
#'
#' @param ... One or more expressions or list of expressions
#'
#' @return A list of expressions
#'
#' @keywords quo
#' @family quo
#'
#' @export
expr_c <- function(...) {
  # Transform single expression into list of expression
  inputs <- map(
    list(...),
    function(x) {
      if (typeof(x) != "list") {
        list(x)
      } else {
        x
      }
    }
  )
  inputs <- flatten(inputs)
  stopifnot(all(map_lgl(inputs, is_expression)))
  is_null <- map_lgl(inputs, is.null)
  inputs[!is_null]
}

#' Replace Expression Value with Name
#'
#' @param expressions A list of expressions
#'
#' @param quosures *Deprecated*, please use `expressions` instead.
#'
#'
#' @keywords quo
#' @family quo
#'
#'
#' @return A list of expressions
#' @export
#'
#' @examples
#' library(rlang)
#' replace_values_by_names(exprs(AVAL, ADT = convert_dtc_to_dt(EXSTDTC)))
replace_values_by_names <- function(expressions, quosures) {
  if (!missing(quosures)) {
    deprecate_stop(
      "0.3.0",
      "replace_values_by_names(quosures = )",
      "replace_values_by_names(expressions = )"
    )
  }
  if (is.null(names(expressions))) {
    return(expressions)
  }
  map2(expressions, names(expressions), function(e, n) {
    if (n == "") {
      return(e)
    }
    as.symbol(n)
  })
}

#' Replace Symbols in an Expression
#'
#' Replace symbols in an expression
#'
#' @param expression Expression
#'
#' @param target Target symbol
#'
#' @param replace Replacing symbol
#'
#' @author Stefan Bundfuss
#'
#' @return The expression where every occurrence of the symbol `target` is
#'   replaced by `replace`
#'
#' @keywords quo
#' @family quo
#'
#' @export
#'
#' @examples
#'
#' library(rlang)
#'
#' replace_symbol_in_expr(expr(AVAL), target = AVAL, replace = AVAL.join)
#' replace_symbol_in_expr(expr(AVALC), target = AVAL, replace = AVAL.join)
#' replace_symbol_in_expr(expr(desc(AVAL)), target = AVAL, replace = AVAL.join)
replace_symbol_in_expr <- function(expression,
                                   target,
                                   replace) {
  assert_expr(expression)
  target <- assert_symbol(enexpr(target))
  replace <- assert_symbol(enexpr(replace))
  if (is.symbol(expression)) {
    if (expression == target) {
      expression <- replace
    }
  } else {
    for (i in seq_along(expression)) {
      if (expression[[i]] == target) {
        expression[[i]] <- replace
      }
    }
  }
  expression
}

#' Add a Suffix to Variables in a List of Expressions
#'
#' Add a suffix to variables in a list of expressions
#'
#' @param order List of expressions
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'   created by `exprs()`, e.g., `exprs(ADT, desc(AVAL))`
#'
#' @param vars Variables to change
#'
#'   *Permitted Values*: list of variables created by `exprs()`
#'
#' @param suffix Suffix
#'
#'   *Permitted Values*: A character scalar
#'
#'
#' @return The list of expression where for each element the suffix (`suffix`) is
#'   added to every symbol specified for `vars`
#'
#' @keywords quo
#' @family quo
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#'
#' add_suffix_to_vars(exprs(ADT, desc(AVAL), AVALC), vars = exprs(AVAL), suffix = ".join")
add_suffix_to_vars <- function(order,
                               vars,
                               suffix) {
  assert_expr_list(order)
  assert_vars(vars)
  assert_character_scalar(suffix)
  for (i in seq_along(vars)) {
    order <- lapply(
      order,
      replace_symbol_in_expr,
      target = !!vars[[i]],
      replace = !!sym(paste0(as_label(vars[[i]]), suffix))
    )
  }
  order
}
