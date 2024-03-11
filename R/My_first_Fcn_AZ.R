#' @name my_first_fcn
#' @param hw true or false
#' @title Prints a welcome message to the Admiral family.
#' @return character vector containing the welcome message.
#' @export

my_first_fcn <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family!")
  }
}
