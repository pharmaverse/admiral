#' Welcome Message Function
#'
#' Displays a friendly welcome message, specifically designed for the "admiral family".
#'
#' @param hw Logical. Determines whether to display a slightly more enthusiastic welcome.
#'           Defaults to TRUE.
#'
#' @details This is my first trial to write a function
#'
#' @return  None. This function's primary purpose is to print a message to the console.
#'
#' @export
#'
#' @examples
#' my_first_fcn()
#' my_first_fcn(hw = FALSE)
#'
#'
my_first_fcn <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family!")
  }
}

