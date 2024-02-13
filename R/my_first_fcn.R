#' Derive on boarding dummy function
#'
#' @param    first_time TRUE or FALSE
#' @details  outputs "Welcome to the admiral family!".
#' @return   a Welcome message
#' @export
#' @keywords utils_print
#' @family   utils_print
#' @examples
#' my_first_fcn(first_time = FALSE)
my_first_fcn <- function(first_time = TRUE) {
  # Check that first_time is Logical
  stopifnot(is.logical(first_time))

  if (first_time) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome again to the admiral family!")
  }
}
