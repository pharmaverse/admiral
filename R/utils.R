#' @export
rlang::expr

#' @export
rlang::exprs

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
#' Negated Value Matching
#'
#' Returns a `logical` vector indicating if there is *no* match of the
#' left operand in the right operand.
#'
#' @param x The values to be matched
#' @param table The values to be matched against
#'
#' @rdname utils
#' @name not_in
#'
#' @export
#'
#' @examples
#' "a" %!in% c("b", "v", "k")
`%!in%` <- function(x, table) {
  !(x %in% table)
}
