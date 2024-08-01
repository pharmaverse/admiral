#' @name my_first_fcn.R
#'
#' @title Hello Admiral Function
#'
#' @param hw TRUE or FALSE
#'
#' @details First function to generate welcome message
#'
#' @return Happy Message
#'
#' @family der_adxx
#'
#' @keywords der_adxx
#'
#' @export
#'
#' @examples
#' hello_admiral(hw = FALSE)
hello_admiral <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family!")
  }
}

