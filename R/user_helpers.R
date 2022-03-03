#' Open an ADaM Template Script
#'
#' Uses list_all_templates() to discover which templates are available
#'
#' @param adam_name An ADaM dataset name.
#' @param save_path Path to save the script.
#' @param overwrite Whether to overwrite an existing file named `save_path`.
#' @param open Whether to open the script right away.
#'
#' @return No return values, called for side effects
#'
#' @author Shimeng Huang, Thomas Neitmann
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
                            overwrite = FALSE,
                            open = interactive()) {
  assert_character_scalar(adam_name)
  assert_character_scalar(save_path)
  assert_logical_scalar(overwrite)
  assert_logical_scalar(open)

  if (!toupper(adam_name) %in% list_all_templates()) {
    err_msg <- paste0(
      sprintf("No template for '%s' available.\n", toupper(adam_name)),
      "\u2139 Run `list_all_templates()` to get a list of all available ADaM templates."
    )
    abort(err_msg)
  }

  if (file.exists(save_path) && !overwrite) {
    err_msg <- paste(
      sprintf("A file named '%s' already exists.", save_path),
      "\u2139 Set `overwrite = TRUE` to force overwriting it.",
      sep = "\n"
    )
    abort(err_msg)
  }

  template_file <- system.file(
    paste0("templates/ad_", tolower(adam_name), ".R"),
    package = "admiral"
  )

  if (file.copy(template_file, save_path, overwrite = TRUE)) {
    inform(sprintf("\u2713 File '%s' has been created successfully", save_path))
  }

  if (open) {
    utils::file.edit(save_path)
  }

  invisible(TRUE)
}

#' List All Available ADaM Templates
#'
#' @author Shimeng Huang, Thomas Neitmann
#'
#' @keywords user_utility
#'
#' @return A `character` vector of all available templates
#'
#' @export
#'
#' @examples
#' list_all_templates()
list_all_templates <- function() {
  list.files(system.file("templates", package = "admiral")) %>%
    str_remove(".R$") %>%
    str_remove("^ad_") %>%
    toupper() %>%
    structure(class = c("adam_templates", "character"))
}

print.adam_templates <- function(x, ...) {
  cat("Existing templates:\n")
  cat(paste0("\U2022 ", x), sep = "\n")
}
