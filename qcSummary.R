# nolint start
library(diffdf)
library(tibble)
library(stringr)
# base or input: is older ds;
# compare:  is new ds (created by template)
# nolint end
input_path <- file.path("inst/verify/old")
comp_path <- file.path("inst/verify/new")

input_dataset_paths <- list.files(input_path)
input_dataset_paths <- input_dataset_paths[endsWith(input_dataset_paths, ".rds")]
input_dataset_names <- tools::file_path_sans_ext(input_dataset_paths)
comp_dataset_paths <- list.files(comp_path)
comp_dataset_names <- tools::file_path_sans_ext(comp_dataset_paths)
#
## restores files to environment
for (i in seq_along(input_dataset_names)) {
  assign(
    paste0("new_", input_dataset_names[i]),
    readRDS(file.path(input_path, input_dataset_paths[i]))
  )
}
#
for (i in seq_along(comp_dataset_names)) {
  if (!file.exists(file.path(comp_path, comp_dataset_paths[i]))) {
    next
  }
  assign(
    paste0("comp_", comp_dataset_names[i]),
    readRDS(file.path(comp_path, comp_dataset_paths[i]))
  )
}

sink("result.Rmd")
cat("## Verify Templates Check Complete!", "\n\n")
cat("Date: ", format(Sys.Date()), "\n")
cat("Run by: ", Sys.getenv("GITHUB_ACTOR"), "\n")
cat("Git Ref: ", Sys.getenv("GITHUB_REF"), "\n")
cat("BASE: ", "Generated ADaM Datasets from Templates during Run", "\n")
cat("COMPARE: ", "ADaM Datasets from pharmaverseadam ", "\n")
for (y in input_dataset_names) {

  new_dataset <- paste0("new_", y)
  comp_dataset <- paste0("comp_", y)
  diffs <- diffdf(base = get(comp_dataset), compare = get(new_dataset))
  # nolint start
  #  if  (diffdf::diffdf_has_issues(diffs)) print(diffs)
  #  if (length(diffs) != 0) file.create("qc.fail")
  # nolint end
  cat("<details>\n")
  status_emoji <- if (length(diffs) == 0) "✅" else "❌"
  cat(str_glue("<summary>{status_emoji} Dataset: {y}</summary>\n\n"))
  cat("\n\n```\n\n")
  print(diffs)
  cat("```\n\n")
  cat("</details>")
  cat("\n\n")
}
sink()

cat("--- END ---\n")
readLines("result.Rmd") |> cat(sep = "\n")
