#' Derive Extension Example
#'
#' Says hello admiral
#'
#'  @param hw TRUE or FALSE
#'
#'  @author author sc460144
#'
#' @details welcome to admiral family message
#'
#' @return Happy Message
#'
#' @family
#'
#' @keywords
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

