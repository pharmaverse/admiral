#' Dummy function for on boarding
#'
#' @param input - required input information in any format
#' @details create a function to output a text message 'Welcome to the admiral family!'
#' @return a welcome message
#' @keywords der_gen
#' @family der_gen
#' @export
#' @examples
#' my_first_fcn("Y")
#' my_first_fcn(NULL)
my_first_fcn <- function(input) {
  if (!is.null(input)) {
    print("Welcome to the admiral family!")
  } else {
    print("Incorrect input")
  }
}
