#' Open an ADaM Template Script
#'
#' @param adam_name An ADaM dataset name. You can use any of the available
#'   dataset names
#'   `r map_chr(list_all_templates(), ~ paste0("\\code{\"", .x, "\"}"))`.
#'   The dataset name is case-insensitive. The default dataset name is `"ADSL"`.
#' @param save_path Path to save the script.
#' @param package The R package in which to look for templates. By default `"admiral"`.
#' @param overwrite Whether to overwrite an existing file named `save_path`.
#' @param open Whether to open the script right away.
#'
#' @return No return values, called for side effects
#'
#' @details Running without any arguments such as `use_ad_template()` auto-generates `adsl.R` in
#' the current path. Use `list_all_templates()` to discover which templates are available.
#'
#' @family utils_examples
#' @keywords utils_examples
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
    cli_abort(c(
      "No template for {toupper(adam_name)} available in package {.pkg {package}}.",
      i = paste(
        "Run {.run admiral::list_all_templates(\"{package}\")} to get a list of",
        "all available ADaM templates."
      )
    ))
  }

  if (file.exists(save_path) && !overwrite) {
    cli_abort(c(
      "A file named {.file {save_path}} already exists.",
      i = "Set {.code overwrite = TRUE} to force overwriting it."
    ))
  }

  template_file <- system.file(
    paste0("templates/ad_", tolower(adam_name), ".R"),
    package = package
  )

  if (file.copy(template_file, save_path, overwrite = TRUE)) {
    cli_inform(c(v = "File {.file {save_path}} has been created successfully"))
  }

  if (open) {
    file.edit(save_path) # nocov
  }

  invisible(TRUE)
}

#' List All Available ADaM Templates
#'
#' @param package The R package in which to look for templates. By default `"admiral"`.
#'
#'
#' @family utils_examples
#' @keywords utils_examples
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
    cli_abort(
      "No package called {.pkg {package}} is installed and hence no templates are available."
    )
  }

  list.files(system.file("templates", package = package)) %>%
    str_remove(".R$") %>%
    str_remove("^ad_") %>%
    toupper() %>%
    structure(class = c("adam_templates", "character"), package = package) # nolint: undesirable_function_linter
}

#' Print `adam_templates` Objects
#'
#' @param x A `adam_templates` object
#' @param ... Not used
#'
#' @return No return value, called for side effects
#'
#'
#' @export
#'
#' @keywords internal
#' @family utils_print
#'
#' @seealso [list_all_templates()]
#'
#' @examples
#' templates <- list_all_templates()
#' print(templates)
print.adam_templates <- function(x, ...) {
  pkg <- attr(x, "package")
  if (length(x) == 0L) {
    cat("No ADaM templates available in package '", pkg, "'\n", sep = "")
  } else {
    cat("Existing ADaM templates in package '", pkg, "':\n", sep = "")
    bullet <- if (is.na(iconv("\U2022"))) "-" else "\U2022"
    cat(paste(bullet, x), sep = "\n")
  }
}
