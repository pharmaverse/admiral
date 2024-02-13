#' Derive onboarding dummy function
#'
#' @author   Fanny Gautier
#'
# Mandatory fields in header:
#' @param    first_time TRUE or FALSE
#' @details  outputs "Welcome to the admiral family!".
#' @return   a Welcome message
#' @export
#' @keywords utils_print
#' @family   utils_print
#  keywords  der_message is not a known tag, so not included for this function
#  family    der_message is not a known tag, so not included for this function
#' @examples my_first_fcn(first_time = FALSE)

my_first_fcn <- function(first_time = TRUE) {
  # Check that first_time is Logical
  expect(is.logical(first_time), "Error: first_time must be logical")

  if (first_time) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome again to the admiral family!")
  }
}
