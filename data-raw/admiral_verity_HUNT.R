
## PURPOSE:    use diffdf to  find which column in adpc is causing problems
library(diffdf)
library(tibble)
library(stringr)
input_path <- "inst/verify/old"
comp_path <- "inst/verify/new"



# load data   (did you run 1st part?)
input_dataset_paths <- list.files(input_path)
input_dataset_paths <- input_dataset_paths[endsWith(input_dataset_paths, ".rds") ]
input_dataset_names <- tools::file_path_sans_ext(input_dataset_paths)
comp_dataset_paths <- list.files(comp_path)
comp_dataset_names <- tools::file_path_sans_ext(comp_dataset_paths)

# restores files to environment
for (i in seq_along(input_dataset_names)) {
  assign(
    paste0("new_", input_dataset_names[i]), 
    readRDS(file.path(input_path, input_dataset_paths[i]))
    )
}

for (i in seq_along(comp_dataset_names)) {
  if (!file.exists(file.path(comp_path, comp_dataset_paths[i]))) {
    next
  }
  assign(
    paste0("comp_", comp_dataset_names[i]), 
    readRDS(file.path(comp_path, comp_dataset_paths[i]))
    )
}

# compare comp_adpc (original) and new_adpc (just created)
identical(comp_adpc, new_adpc)   # FALSE

old = comp_adpc
new = new_adpc

dim(old)  # 4479 x 127
dim(new)  # 3852 x 127


# Method A, download raw from github
# dim 4479 x 127


# Method B, as now done = WRONG

libary(pharmaverseadam)
adam_names=c("adpc")
path = "inst/verify/old"

  f = function(ds) {
      q = call("::", "pharmaverseadam", sym(ds))
      save(q,
          file = file.path(path, paste0(ds, ".rda"))) }

  walk(.x = adam_names, .f = f)

B1 = get(load("inst/verify/old/adpc.rda"))
eval(B1)   # 3852 x 127


#' Copies ADaM datasets from pharmaverseadam
#' @param adam_names character vector  Set of  ADaMs to download.
#' @param path Character string. Directory to save downloaded ADaMs.
download_adam_old <- function(adam_names, path = NULL) {

# method C, legacy
adam_names="adpc"
adam = "adpc"
path = "inst/verify/old"

   lapply(adam_names, function(adam) {
     githubURL <- paste0( # nolint
       "https://github.com/pharmaverse/pharmaverseadam/raw/refs/heads/main/data/",
       adam, ".rda?raw=true"
     )
     cat("Downloading: ", adam, "\n")
     download.file(
       url = githubURL,
       quiet = TRUE,
       destfile = paste0(path, "/", adam, ".rda"),
       mode = "wb"
     )
   })

dim(get(load("inst/verify/old/adpc.rda")))   # 4479 x 127


# method C
pak::pak("pharmaverseadam")
load_all()
C = pharmaverseadam::adpc
dim(C)   # 3852 x 127

# D cloned, local copy
# D = load("~/code/pharmaverseadam/data/adpc.rda")

key=c("USUBJID", "PARAMCD", "AVISTIN", "ATPTN", "DTYPE")

sink("qc.Rmd")
cat("## Verify Templates Check Complete!", "\n\n")
cat("Date: ", format(Sys.Date()), "\n")
cat("Run by: ", Sys.getenv('GITHUB_ACTOR'), "\n")
cat("Git Ref: ", Sys.getenv("GITHUB_REF"), "\n")
cat("BASE: ", "Generated ADaM Datasets from Templates during Run", "\n")
cat("COMPARE: ", "ADaM Datasets from pharmaverseadam ", "\n")
for (y in input_dataset_names) {
  if (y == "adpc") next()
  new_dataset <- paste0("new_", y)
  comp_dataset <- paste0("comp_", y)
  diffs <- diffdf(base = get(comp_dataset), compare =  get(new_dataset))
#  if  (diffdf::diffdf_has_issues(diffs)) print(diffs)
    
#  if (length(diffs) != 0) file.create("qc.fail")
    
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
```

```{r data_check}
#| warning: false
#| message: false
#| echo: false
#| output: asis
readLines("qc.Rmd") |>
  cat(sep = "\n")
```
