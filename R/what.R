#' What Kind of Object is This?
#'
#' Returns a string describing what kind of object the input is.
#'
#' @param x Any R object
#'
#' @author Thomas Neitmann
#'
#' @keywords what
#' @family what
#'
#' @rdname dev_util_what_is_it
#'
#' @examples
#' admiral:::what_is_it(mtcars)
#' admiral:::what_is_it(NA)
#' admiral:::what_is_it(TRUE)
#' admiral:::what_is_it(lm(hp ~ mpg, data = mtcars))
#' admiral:::what_is_it(letters)
what_is_it <- function(x) {
  if (is.null(x)) {
    "`NULL`"
  } else if (is.factor(x)) {
    "a factor"
  } else if (is.symbol(x)) {
    "a symbol"
  } else if (isS4(x)) {
    sprintf("a S4 object of class '%s'", class(x)[1L])
  } else if (is.atomic(x) && length(x) == 1L) {
    if (is.character(x)) {
      paste0("`\"", x, "\"`")
    } else {
      paste0("`", x, "`")
    }
  } else if (is.atomic(x) || class(x)[1L] == "list") {
    friendly_type_of(x)
  } else if (is.data.frame(x)) {
    "a data frame"
  } else {
    sprintf("an object of class '%s'", class(x)[1L])
  }
}
