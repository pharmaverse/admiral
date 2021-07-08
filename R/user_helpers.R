#' Open an ADaM Template Script
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
  assert_character_scalar(adam_name)
  assert_character_scalar(save_path)
  assert_logical_scalar(open)

  if (!requireNamespace("usethis", quietly = TRUE)) {
    abort("Required package {usethis} is not installed.")
  }

  if (!toupper(adam_name) %in% list_all_templates()) {
    err_msg <- paste0(
      sprintf("No template for '%s' available.\n", toupper(adam_name)),
      "ℹ Run `list_all_templates()` to get a list of all available ADaM templates."
    )
    abort(err_msg)
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
  list.files(system.file("templates", package = "admiral")) %>%
    str_remove(., ".R$") %>%
    str_remove(., "^ad_") %>%
    toupper(.) %>%
    structure(class = c("adam_templates", "character"))
}

print.adam_templates <- function(x, ...) {
  cat("Existing templates:\n")
  cat(paste0("\U2022 ", x), sep = "\n")
}
