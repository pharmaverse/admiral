#' Concatenate One or More Quosure(s)
#'
#' @param ... One or more objects of class `quosure` or `quosures`
#'
#' @return An object of class `quosures`
#'
#' @author Thomas Neitmann
#'
#' @keywords quo
#' @family quo
#'
#' @export
quo_c <- function(...) {
  inputs <- unlist(list(...), recursive = TRUE)
  stopifnot(all(map_lgl(inputs, is_quosure)))
  is_null <- map_lgl(inputs, quo_is_null)
  rlang::as_quosures(inputs[!is_null])
}

#' Check Whether an Argument Is Not a Quosure of a Missing Argument
#'
#' @param x Test object
#'
#' @return TRUE or error.
#'
#' @author Thomas Neitmann, Ondrej Slama
#'
#' @keywords quo
#' @family quo
#'
#' @export
quo_not_missing <- function(x) {
  !rlang::quo_is_missing(x)

  if (is.null(missing(x)) || quo_is_missing(x)) {
    stop(paste0(
      "Argument `",
      deparse(substitute(x)),
      "` is missing, with no default"
    ))
  }
}


#' Replace Quosure Value with Name
#'
#' @param quosures A list of quosures
#'
#' @author Thomas Neitmann
#'
#' @keywords quo
#' @family quo
#'
#'
#' @return A list of quosures
#' @export
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

#' Replace Symbols in a Quosure
#'
#' Replace symbols in a quosure
#'
#' @param quosure Quosure
#'
#' @param target Target symbol
#'
#' @param replace Replacing symbol
#'
#' @author Stefan Bundfuss
#'
#' @return The quosure where every occurence of the symbol `target` is replaced
#'   by `replace`
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
#' replace_symbol_in_quo(quo(AVAL), target = AVAL, replace = AVAL.join)
#' replace_symbol_in_quo(quo(AVALC), target = AVAL, replace = AVAL.join)
#' replace_symbol_in_quo(quo(desc(AVAL)), target = AVAL, replace = AVAL.join)
replace_symbol_in_quo <- function(quosure,
                                  target,
                                  replace) {
  assert_expr(quosure)
  target <- quo_get_expr(assert_symbol(enquo(target)))
  replace <- quo_get_expr(assert_symbol(enquo(replace)))
  expr <- quo_get_expr(quosure)
  if (is.symbol(expr)) {
    if (expr == target) {
      expr <- replace
    }
  } else {
    for (i in seq_along(quosure)) {
      if (expr[[i]] == target) {
        expr[[i]] <- replace
      }
    }
  }
  rlang::quo_set_expr(quosure, expr)
}

#' Add a Suffix to Variables in a List of Quosures
#'
#' Add a suffix to variables in a list of quosures
#'
#' @param order List of quosures
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'   created by `vars()`, e.g., `vars(ADT, desc(AVAL))`
#'
#' @param vars Variables to change
#'
#'   *Permitted Values*: list of variables created by `vars()`
#'
#' @param suffix Suffix
#'
#'   *Permitted Values*: A character scalar
#'
#' @author Stefan Bundfuss
#'
#' @return The list of quosures where for each element the suffix (`suffix`) is
#'   added to every symbol specified for `vars`
#'
#' @keywords quo
#' @family quo
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' add_suffix_to_vars(vars(ADT, desc(AVAL), AVALC), vars = vars(AVAL), suffix = ".join")
add_suffix_to_vars <- function(order,
                               vars,
                               suffix) {
  assert_order_vars(order)
  assert_vars(vars)
  assert_character_scalar(suffix)
  for (i in seq_along(vars)) {
    order <- lapply(
      order,
      replace_symbol_in_quo,
      target = !!quo_get_expr(vars[[i]]),
      replace = !!sym(paste0(as_label(
        quo_get_expr(vars[[i]])
      ), suffix))
    )
  }
  class(order) <- c("quosures", "list")
  order
}
