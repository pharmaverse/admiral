#' Is a named argument
#'
#' @param x Any R object
#'
#' @return `TRUE` if the argument is named, `FALSE` otherwise
#' @export
#'
#' @keywords is
#' @family is
#' @export
is_named <- function(x) {
  !is.null(names(x)) && all(names(x) != "")
}

#' Checks if the argument equals the auto keyword
#'
#' @param arg argument to check
#'
#' @return `TRUE` if the argument equals the auto keyword, i.e., it is a quosure
#'   of a symbol named auto.
#'
#' @author Stefan Bundfuss
#'
#' @keywords is
#' @family is
#' @export
is_auto <- function(arg) {
  is_quosure(arg) && quo_is_symbol(arg) && quo_get_expr(arg) == expr(auto)
}

#' Is order vars?
#'
#' Check if inputs are created using `vars()` or calls involving `desc()`
#' @param arg An R object
#'
#' @return `FALSE` if the argument is not a list of order vars
#'
#' @export
#'
#' @keywords is
#' @family is
is_order_vars <- function(arg) {
  quo_is_desc_call <- function(quo) {
    expr <- quo_get_expr(quo)
    is_call(expr) &&
      length(expr) == 2L &&
      deparse(expr[[1L]]) == "desc" &&
      is_symbol(expr[[2L]])
  }

  inherits(arg, "quosures") &&
    all(map_lgl(arg, ~ quo_is_symbol(.x) || quo_is_desc_call(.x)))
}

#' Is this string a valid DTC
#'
#' @param arg A `character` vector
#'
#' @return `TRUE` if the argument is a valid `--DTC` string, `FALSE` otherwise
#' @export
#' @keywords is
#' @family is
#'
is_valid_dtc <- function(arg) {
  twod <- "(\\d{2}|-)"
  pattern <- paste0(
    "^(\\d{4}|-)?",
    "(-", twod, ")?",
    "(-", twod, ")?",
    "(T", twod, ")?",
    "(:", twod, ")?",
    "(:", twod, "(.(\\d{1,5}))?)?$"
  )

  grepl(pattern, arg) | arg == "" | is.na(arg)
}
