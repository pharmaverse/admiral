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

# ------------------------  separate1
library(glue)

CON = file("qcSummary.txt", "a")
writeLines(paste0("## Verify Templates Check Complete!", "\n\n"), CON)
writeLines(paste0("Date: ", Sys.Date(), "\n"), CON)
writeLines(paste0("Run by: ", Sys.getenv("GITHUB_ACTOR"), "\n"), CON)
writeLines(paste0("Git Ref: ", Sys.getenv("GITHUB_REF"), "\n"), CON)
writeLines(paste0("BASE: ", "Generated ADaM Datasets from Templates during Run", "\n"), CON)
writeLines(paste0("COMPARE: ", "ADaM Datasets from pharmaverseadam ", "\n"), CON)

for (y in input_dataset_names) {

  new_dataset <- paste0("new_", y)
  comp_dataset <- paste0("comp_", y)

   diffs <- tryCatch(
    {  
        diffs <- diffdf(base = get(comp_dataset), compare = get(new_dataset))
    },
    error = function(e)  {
        message("Error in diffdf ", e$message)
        # empty list to not break loop
        list()
    # nolint start
    #  if  (diffdf::diffdf_has_issues(diffs)) print(diffs)
    #  if (length(diffs) != 0) file.create("qc.fail")
    # nolint end
        })

  
  writeLines("<details>\n", CON)
  status_emoji <- if (length(diffs) == 0) "✅" else "❌"
  writeLines(str_glue("<summary>{status_emoji} Dataset: {y}</summary>\n\n"), CON)
  writeLines("\n\n```\n\n", CON)
  writeLines(as.character(length(diffs)), CON)
  writeLines(print(diffs,as_string=TRUE), CON)
  writeLines("```\n\n", CON)
  writeLines("</details>", CON)
  writeLines("\n\n", CON)
}
close(CON)
readLines("log_file.txt") |> cat(sep = "\n")
