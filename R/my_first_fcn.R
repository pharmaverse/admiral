#' Welcome message
#' @param hw TRUE or FALSE
#' @return A message to the user
#' @export
#'
#' @examples
#' my_first_fcn()
my_first_fcn <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family")
  }
}
