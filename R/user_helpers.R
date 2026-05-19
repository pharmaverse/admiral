#' Open an ADaM Template Script
#'
#' @param adam_name An ADaM dataset name. You can use any of the available
#'   dataset names
#'   `r map_chr(list_all_templates(), ~ paste0("\\code{\"", .x, "\"}"))`.
#'   The dataset name is case-insensitive. The default dataset name is `"ADSL"`.
#' @param save_path Path to save the script.
#' @param package The R package in which to look for templates. By default `"admiral"`.
#' @param overwrite Whether to overwrite an existing file named `save_path`.
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
                            overwrite = FALSE) {
  assert_character_scalar(adam_name)
  assert_character_scalar(save_path)
  assert_character_scalar(package)
  assert_logical_scalar(overwrite)

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

#' Output Boilerplate for the "Add Labels and Attributes" Vignette Section
#'
#' Outputs the standard markdown text used in ADaM-specific vignettes for the
#' "Add Labels and Attributes" section. This function is intended to be called
#' inside an R Markdown code chunk with `results='asis'` and `echo=FALSE`.
#'
#' @return No return value. The function outputs text directly to the console
#'   or output stream via `cat()`, intended for use in R Markdown documents
#'   with `results='asis'`.
#'
#' @details
#' The outputted section describes how to add variable labels and other
#' metadata to ADaM datasets as a final step in the derivation process,
#' using the `{metacore}`, `{metatools}`, and `{xportr}` packages.
#'
#' @keywords internal
#'
#' @examples
#' admiral_labels_attrs_section()
admiral_labels_attrs_section <- function() {
  cat(paste0(
    "## Add Labels and Attributes {#attributes}\n\n",
    "Note that attributes may not be preserved in some cases after processing\n",
    "with `{admiral}`. The recommended approach is to apply variable labels\n",
    "and other metadata as a final step in your data derivation process using\n",
    "packages like:\n\n",
    "-   [metacore](https://atorus-research.github.io/metacore/): establish a\n",
    "    common foundation for the use of metadata within an R session.\n\n",
    "-   [metatools](https://pharmaverse.github.io/metatools/): enable the\n",
    "    use of metacore objects. Metatools can be used to build datasets or\n",
    "    enhance columns in existing datasets as well as checking datasets\n",
    "    against the metadata.\n\n",
    "-   [xportr](https://atorus-research.github.io/xportr/): functionality\n",
    "    to associate all metadata information to a local R data frame,\n",
    "    perform data set level validation checks and convert into a\n",
    "    [transport v5 file(xpt)](https://documentation.sas.com/doc/en/",
    "pgmsascdc/9.4_3.5/movefile/n1xbwdre0giahfn11c99yjkpi2yb.htm).\n\n",
    "NOTE: Together with `{admiral}` these packages comprise an End to End\n",
    "pipeline under the umbrella of the\n",
    "[pharmaverse](https://github.com/pharmaverse). An example of applying\n",
    "metadata and performing associated checks can be found at the [pharmaverse\n",
    "E2E example](https://pharmaverse.github.io/examples/adam/adsl).\n"
  ))
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
