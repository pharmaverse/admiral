#' A helper function to create a template script
#'
#' @param adam_name An ADaM dataset name.
#'
#' @param save_path Path to save the script.
#'
#' @param open Whether to open the script right away.
#'
#' @author Shimeng Huang
#'
#' @export
#'
#' @examples
#' \dontrun{
#' use_ad_template("./ad_adsl.R", "adsl")
#' }
use_ad_template <- function(adam_name = "adsl",
                            save_path = paste0("./", adam_name, ".R"),
                            open = interactive()) {

  if (!requireNamespace("usethis", quiet = TRUE)) {
    abort("Required package {usethis} is not installed.")
  }

  usethis::use_template(
    template = paste0("ad_", tolower(adam_name), ".R"),
    save_as = save_path,
    package = "admiral",
    open = open
  )
}

#' List all templates provided by {admiral}
#'
#' @author Shimeng Huang
#'
#' @export
#'
#' @examples
#' list_all_templates()
list_all_templates <- function() {
  all_tpl <- list.files(system.file("templates", package = "admiral"))
  all_tpl <- str_remove(all_tpl, ".R$")
  all_tpl <- str_remove(all_tpl, "^ad_")
  all_tpl <- paste0("\U2022 ", toupper(all_tpl))
  cat("Existing templates: \n")
  cat(all_tpl, sep = "\n")
}
