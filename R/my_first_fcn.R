#' Derive Extension Example
#'
#' Says hello admiral
#'
#' @param hw TRUE or FALSE
#'
#' @author Anders
#' @keywords other_advanced
#'
#' @return Happy Message
#'
#' @export
#'
#' @examples
#' hello_admiral(hw = FALSE)
hello_admiral <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family hw!")
  } else {
    message("Welcome to the admiral family!")
  }
}
