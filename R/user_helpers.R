#' Open an ADaM Template Script
#'
#' @param adam_name An ADaM dataset name. You can use any of the available dataset name `r list_all_templates()`, and the dataset name is case-insensitive. The default dataset name is ADSL.
#' @param save_path Path to save the script.
#' @param overwrite Whether to overwrite an existing file named `save_path`.
#' @param open Whether to open the script right away.
#'
#' @return No return values, called for side effects
#'
#' @details Running without any arguments such as `use_ad_template()` auto-generates adsl.R in the current path. Use `list_all_templates()` to discover which templates are available.
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
                            package = "admiral",
                            overwrite = FALSE,
                            open = interactive()) {
  assert_character_scalar(adam_name)
  assert_character_scalar(save_path)
  assert_character_scalar(package)
  assert_logical_scalar(overwrite)
  assert_logical_scalar(open)

  if (!toupper(adam_name) %in% list_all_templates(package)) {
    err_msg <- sprintf(
      "No template for '%s' available in package '%s'.\n",
      "\u2139 Run `list_all_templates('%s')` to get a list of all available ADaM templates.",
      toupper(adam_name), package, package
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
    package = package
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
list_all_templates <- function(package = "admiral") {
  assert_character_scalar(package)

  if (!requireNamespace(package, quietly = TRUE)) {
    err_msg <- sprintf("No package called '%s' is installed and hence no templates are available", package)
    abort(err_msg)
  }

  list.files(system.file("templates", package = package)) %>%
    str_remove(".R$") %>%
    str_remove("^ad_") %>%
    toupper() %>%
    structure(class = c("adam_templates", "character"))
}

print.adam_templates <- function(x, ...) {
  cat("Existing templates:\n")
  cat(paste0("\U2022 ", x), sep = "\n")
}
