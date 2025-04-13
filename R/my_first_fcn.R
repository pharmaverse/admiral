#' Display a Welcome Message to the Admiral Family
#'
#' @description
#' This function displays a welcome message to users of the admiral package.
#' It is for dummy testing.
#'
#' @param hw
#'   `TRUE` to include "hw" in the welcome message, `FALSE` otherwise.
#'   Defaults to `TRUE`.
#'
#' @return
#'   Happy message
#'
#' @family der_adxx
#' @keywords der_adxx
#' @export
#'
#' @examples
#' # Display welcome message with hw
#' hello_admiral()
#'
#' # Display welcome message without hw
#' hello_admiral(hw = FALSE)
#'
hello_admiral <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family with hw")
  } else {
    message("Welcome to the admiral family without hw")
  }
}
