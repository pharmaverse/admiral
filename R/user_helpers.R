#' A helper function to create a template script.
#'
#' @param save_path Path to save the script.
#' @param dom_name An ADaM domain name.
#' @param open Whether to open the script right away.
#'
#' @author Shimeng Huang
#'
#' @importFrom usethis use_template
#'
#' @export
#'
#' @examples
#' use_ad_template(".", "adsl")
use_ad_template <- function(save_path, dom_name = "adsl", open = interactive()) {
  usethis::use_template(
    template = paste0("ad_", tolower(dom_name), ".R"),
    save_as = save_path,
    package = "admiral",
    open = open
  )
}
