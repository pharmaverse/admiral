#' @export
rlang::expr

#'@export
rlang::exprs

#' Enumerate Multiple Strings
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' enumerate(letters[1:6])
enumerate <- function(x) {
  paste(
    paste0(backquote(x[-length(x)]), collapse = ", "),
    "and",
    backquote(x[length(x)])
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
