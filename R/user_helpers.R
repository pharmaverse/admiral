#' Open a ADaM Template Script
#'
#' @param adam_name An ADaM dataset name.
#'
#' @param save_path Path to save the script.
#'
#' @param open Whether to open the script right away.
#'
#' @author Shimeng Huang
#'
#' @keywords user_utility
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'   use_ad_template("adsl")
#' }
use_ad_template <- function(adam_name = "adsl",
                            save_path = paste0("./", adam_name, ".R"),
                            open = interactive()) {
  if (!requireNamespace("usethis", quietly = TRUE)) {
    abort("Required package {usethis} is not installed.")
  }

  usethis::use_template(
    template = paste0("ad_", tolower(adam_name), ".R"),
    save_as = save_path,
    package = "admiral",
    open = open
  )
}

#' List All Available ADaM Templates
#'
#' @author Shimeng Huang
#'
#' @keywords user_utility
#'
#' @export
#'
#' @examples
#' list_all_templates()
list_all_templates <- function() {
  all_tpl <- list.files(system.file("templates", package = "admiral")) %>%
    str_remove(., ".R$") %>%
    str_remove(., "^ad_") %>%
    toupper(.) %>%
    paste0("\U2022 ", .)
  cat("Existing templates: \n")
  cat(all_tpl, sep = "\n")
}
