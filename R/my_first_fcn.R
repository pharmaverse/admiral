#' My first function, happy greeting
#'
#' @param hw TRUE or FALSE
#'
#' @details In the roxygen documentation you will find tags for family and keywords.
#' This is to create organized sections for the Reference tab on the pkgdown website.
#' You can modify the `_pkgdown.yml` as necessary to create appropriate sections as necessary.
#' Under `./man/roxygen/meta.R`, you will find where to store these family/keywords.
#'
#' @return Happy Message
#'
#' @family example
#'
#' @keywords example
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
