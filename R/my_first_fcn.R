#' My first function, happy greeting
#'
#' @param smile TRUE or FALSE
#'
#' @details This is a message of happiness.
#'
#' @return Happy Message
#'
#' @family example
#'
#' @keywords example
#'
#' @export
#'
#' @examples
#' hello_admiral(smile = FALSE)
hello_admiral <- function(smile = TRUE) {
  if (smile) {
    message("Welcome to the admiral family :)")
  } else {
    message("Welcome to the admiral family!")
  }
}
