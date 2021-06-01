#' @export
rlang::expr

#' @export
rlang::exprs

#' @export
dplyr::vars

#' Enumerate Multiple Strings
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' enumerate(letters[1:6])
enumerate <- function(x, quote_fun = backquote) {
  paste(
    paste0(quote_fun(x[-length(x)]), collapse = ", "),
    "and",
    quote_fun(x[length(x)])
  )
}

#' Wrap a String in Backquotes
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' backquote("foo")
backquote <- function(x) {
  paste0("`", x, "`")
}

#' Wrap a String in SIngle Quotes
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' squote("foo")
squote <- function(x) {
  paste0("'", x, "'")
}

`%!in%` <- Negate(`%in%`)
