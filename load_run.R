# standalone R script
# source("load_run.R")

library(diffdf)
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @return
##' @author jim
base <- readRDS("old.RDS")
compare <- readRDS("new.RDS")
keys <- NULL
file <- "diff_output.txt"

res <- diffdf::diffdf(
  base = base,
  compare = compare,
  file = file,
  suppress_warnings = TRUE # for now
)
res

sink(file = "qc.Rmd")
cat("## Dataset QC Check Complete!", "\n\n")
cat("Date: ", Sys.Date(), "\n")
cat("Run by: ", Sys.getenv("GITHUB_ACTOR"), "\n")
cat("Git Ref: ", Sys.getenv("GITHUB_REF"), "\n")
# for (y in input_dataset_names) {
#  new_dataset <- paste0("new_", y)
#  comp_dataset <- paste0("comp_", y)
#  diffs <- diffdf(get(new_dataset), get(comp_dataset))
#
#  cat("<details>\n")
#  status_emoji <- if (length(diffs) == 0) "✅" else "❌"
#  cat(str_glue("<summary>{status_emoji} Dataset: {y}</summary>\n\n"))
#  cat("\n\n```\n\n")
#  print(diffs)
#  cat("```\n\n")
#  cat("</details>")
#   cat("\n\n")
# }
diffdf::diffdf(
  base = readRDS("old.RDS"),
  compare = readRDS("new.RDS"),
  file = "output_diff.txt"
)
print(diffdf)
sink()

#| warning: false
#| message: false
#| echo: false
#| output: asis
readLines("qc.Rmd") |>
  cat(sep = "\n")

sink("jim.txt")
i <- 1:3
print(i)
sink()
