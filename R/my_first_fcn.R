#' Derive Extension Example
#'
#' Says hello admiral
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
