library(purrr)
library(cli)
library(stringr)
library(devtools)
devtools::install_github("pharmaverse/pharmaversesdtm", ref = "main")
devtools::install_github("pharmaverse/pharmaverseadam", ref = "main")
library(pharmaverseadam)
#nolint start

# Gather templates ----
template_dir <- file.path("inst/templates")
templates <- list.files(template_dir, pattern = "ad_")

# From templates generate vector of adam_names ----
adam_names <- vapply(
  templates, function(x) gsub("ad_|\\.R", "", x),
  USE.NAMES = FALSE, character(length = 1)
)

# Ignore ADLBHY
adam_names <- adam_names[adam_names != "adlbhy"]

# Run templates and compare ----
diffs <- purrr::map(adam_names, function(adam) {
  cli_inform("Running {adam} template and comparing with pharmaverseadam version...")

  ## Get {pharmaverseadam} dataset ----
  data_env <- new.env()
  data(list = adam, package = "pharmaverseadam", envir = data_env)
  dataset_old <- data_env[[adam]]

  ## Produce {admiral} dataset ----
  source(paste0("inst/templates/ad_", adam, ".R"), echo = FALSE)
  dataset_new <- get(adam)

  ## Remove column attributes from both ----
  for (name in names(dataset_new)) {
    attr(dataset_new[[name]], "label") <- NULL
  }

  for (name in names(dataset_old)) {
    attr(dataset_old[[name]], "label") <- NULL
  }

  ## Do comparison ----
  res <- diffdf::diffdf(dataset_new, dataset_old, suppress_warnings = TRUE)

  if (diffdf::diffdf_has_issues(res)) {
    cli_inform("Differences identified in {adam}! See diffdf output at the end.")
    res
  } else {
    cli_inform("No differences identified in {adam}.")
  }
}) %>%
  purrr::set_names(adam_names)

cli_inform("Done with running ADaM templates.\n")

# Return differences, if any ----
diffs_with_differences <- diffs[!sapply(diffs, is.null)]

if (length(diffs_with_differences) > 0) {
  print(diffs_with_differences)
  cli::cli_abort("Erroring due differences between admiral and pharmaverseadam templates.")
}

#nolint end
