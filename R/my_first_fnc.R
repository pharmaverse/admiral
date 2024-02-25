#' Greetings
#'
#' @param hw logical
#'
#' @return Happy message
#'
#' @family internal
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' hello_admiral(hw = FALSE)
hello_admiral <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("It's okay, still welcome to the admiral family :)")
  }
}
