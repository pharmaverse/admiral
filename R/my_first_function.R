#' @title hello_admiral
#'
#' @param hw logical (TRUE or FALSE)
#'
#' @details Use roxygen2
#' Under ./man/roxygen/meta.R

#' @return string
#' @family der_adx
#' @keywords der_adxx
#' @examples
#' hello_admiral(hw = TRUE)
#' @export

hello_admiral = function (hw = TRUE) {
    if (hw) {
        message("Welcome to Admiral family")
    } else {
        cli_text("Using cli:: rather than base:: to say   'Welcome to Admiral Family'")
        #message("Welcome to Admiral family")
    }
}
