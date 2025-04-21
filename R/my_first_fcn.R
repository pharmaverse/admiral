#' Welcome message function
#'
#' @description
#' This function displays a welcome message to users of the admiral package and is used for on-boarding purposes.
#'
#' @param exclamation `TRUE` or `FALSE`
#'
#' @details
#' The `exclamation` argument refers to the presence or absence of an exclamation point (e.g., "!") in the output.
#'
#' @return
#'   A string that welcomes the user to the admiral family, with or without an exclamation point.
#'
#' @family utils_examples
#'
#' @keywords utils_examples
#'
#' @export
#'
#' @examples
#' # Display welcome message without an exclamation point
#' hello_admiral(exclamation = FALSE)
#'
#' # Display welcome message with an exclamation point
#' hello_admiral()
#'
hello_admiral <- function(exclamation = TRUE) {
  if (exclamation) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family.")
  }
}
