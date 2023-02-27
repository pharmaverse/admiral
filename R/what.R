#' What Kind of Object is This?
#'
#' Returns a string describing what kind of object the input is.
#'
#' @param x Any R object
#'
#' @return A `character` description of the type of `x`
#'
#'
#' @keywords what
#' @family what
#'
#' @export
#'
#' @examples
#' what_is_it("abc")
#' what_is_it(1L)
#' what_is_it(1:10)
#' what_is_it(mtcars)
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
