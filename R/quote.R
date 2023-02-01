#' Enumerate Multiple Elements
#'
#' Enumerate multiple elements of a vector or list.
#'
#' @param x A vector or list
#' @param quote_fun Quoting function, defaults to `backquote`. If set to `NULL`,
#'   the elements are not quoted.
#' @param conjunction Character to be used in the message, defaults to `"and"`.
#'
#'
#' @return A `character` vector
#'
#' @keywords quote
#' @family quote
#'
#' @examples
#' enumerate(c("one", "two", "three"))
#'
#' enumerate(c(1, 2, 3), quote_fun = NULL)
#'
#' @export
enumerate <- function(x, quote_fun = backquote, conjunction = "and") {
  if (is.null(quote_fun)) {
    quote_fun <- function(x) x
  }
  if (length(x) == 1L) {
    quote_fun(x)
  } else {
    paste(
      paste0(quote_fun(x[-length(x)]), collapse = ", "),
      conjunction,
      quote_fun(x[length(x)])
    )
  }
}

#' Wrap a String in Backquotes
#'
#' @param x A `character` vector
#'
#'
#' @return A `character` vector
#'
#' @keywords quote
#' @family quote
#'
#' @export
backquote <- function(x) {
  paste0("`", x, "`")
}

#' Wrap a String in Single Quotes
#'
#' @param x A `character` vector
#'
#'
#' @return A `character` vector
#'
#' @keywords quote
#' @family quote
#'
#' @export
squote <- function(x) {
  paste0("'", x, "'")
}

#' Wrap a String in Double Quotes
#'
#' Wrap a string in double quotes, e.g., for displaying character values in
#' messages.
#'
#' @param x A character vector
#'
#' @return If the input is `NULL`, the text `"NULL"` is returned. Otherwise, the
#'   input in double quotes is returned.
#'
#'
#' @keywords quote
#' @family quote
#'
#' @export
dquote <- function(x) {
  if (is.null(x)) {
    "NULL"
  } else {
    paste0("\"", x, "\"")
  }
}
