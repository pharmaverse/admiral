#' Calls a Function Provided by the User
#'
#' Calls a function provided by the user and adds the function call to the error
#' message if the call fails.
#'
#' @param call Call to be executed
#'
#'
#' @return The return value of the function call
#'
#' @family utils_help
#' @keywords utils_help
#'
#' @export
#'
#' @examples
#' hello_admiral(hw=FALSE)
 hello_admiral <- function(hw = TRUE) {
   if (hw) {
     message("My first function in Admiral.")
     } else {
     message("My first function in Admiral.")
     }
}
