# This script:  data-raw/admiral_verify_templates.R

# Assumptions/Questions:
# - ignore *.rda files in admiral/data (per Ben)
# - compares full ADaM - ie all rows
# - for developer use and developer has run load_all()
# - use cli:: for messages/errors - YES
# - most datasets are  stored as .rda files


#' Directories used to find files:

#'  cache_dir
#'  diff_dir
#'  tempdir()
#'  template_dir  location of template R files (inst/templates)
#'  adam_new_dir (aka tempdir() or cache_dir ):  after running, templates place new ADaMs here
#'  adam_old_dir : ADaMs downloaded from pharamverseadam

#' Named vectors/lists, use adam_name to name
#'
#'  adam_names (named character vector)
#'  path  (named list)
#'  keys  (named list)
#'  obj   (named list)
#'  templates (named list)


#' @param adams_names character vector. ADaM or CDISC name, without prefix or suffix  (ex:  adlb)

## DISCUSS, several named vectors, use single tibble as lookup?
#' @param template_dir  Path to templates in the active package (load_all())
#'   ("inst/directory)
#' @param cache_dir (cache)  Path to cache.       (varies by OS/configuration)
#'   ~/.config/R/admiral_templates_data
#' @param adam_old_dir  temporary directory where old ADaMs
#'   (downloaded from github.com) are stored
#' @param adam_new_dir  temporary directory where new ADaMs
#'   (created by templates) are stored


#' @title verify_templates
#' @param pkg package (currently only admiral)
#' @param ds  character vector of ADaM names.  Corresponds to templates to run.  \
#'   A subset of adam_names.
#' @description:
#' Generates ADaM from templates and compares to previously generated ADaM file.
#' (The latter are found in https:://github.com main branch pharmaverseadam).
#'
#' Much code taken from pharamavreseadam::create_adams_data.R
#'  (https://github.com/pharmaverse/pharmaverseadam/blob/main/data-raw/create_adams_data.R)
#'
#' USAGE:   verify_templates()
#'
#' @export
verify_templates <- function(pkg = "admiral", ds = NULL) {
    # ASSUME:  not for interactive use, GHA workflow only

  # no value for ds?  then run all templates
  if (is.null(ds)) {
  ds = c("adae", "adcm", "adeg", "adex", "adlb", "adlbhy", "admh", "adpc",  
           "adpp", "adppk", "adsl", "advs")
  }

  #ds = c("adpp", "adppk", "adsl", "advs")

  pkg <- "admiral"
  if (pkg != "admiral") cli_abort("Currently only `admiral` package is supported")

  # nolint start
  library(pkg, character.only = TRUE)
  library(purrr)
  library(cli)
  library(stringr)
  # nolint end

  cli_inform("Clean directories")


  # always remove, then re-create permanent directories (for .rds)
    if (dir.exists("inst/verify")) {
        res = unlink("inst/verify", recursive=TRUE)}
  dir.create("inst/verify/old", recursive = TRUE)
  dir.create("inst/verify/new", recursive = TRUE)

  # temporary directories:   two cases
  # clean
  if (dir.exists(file.path(tempdir(), "old"))) {
        unlink("*.rda", recursive = TRUE)
        }

  # create fresh (skips if dir  exists)
  path <- create_directories()
            

  cli_alert("Generating ADaMs for { pkg} package.")

  # gather templates ----
  template_dir = file.path("inst/templates")
  templates <- list.files(template_dir, pattern = "ad_")

  cache_dir = tools::R_user_dir("admiral_templates_data", which = "cache")

  # from templates generate vector of adam_names
  adam_names <- vapply(templates, function(x) gsub("ad_|\\.R", "", x),
    USE.NAMES = FALSE, character(length = 1)
  )

  # check
  if (length(templates) != length(adam_names)) cli_abort("Number of templates and adam_names differ")

  # templates is a named chr[]
  names(templates) <- adam_names

  # User specified templates in `ds`. Limit to these templates.
  adam_names <- adam_names[adam_names %in% ds]

  # per Ben, ignore "adlbhy"
  adam_names <- adam_names[adam_names != "adlbhy"]

  # download, saved prior ADaMs from pharmaverseadam as .rda files
  cli_inform("---- Begin copying from github pharmaverseadam")
  download_adam_old(adam_names, path = path$adam_old_dir)

  cli_inform("---- Run templates\n")

    
    purrr::map(adam_names, .progress = TRUE, function(adam){
        cli_inform("------Template running for {adam}")
        ##tryCatch({
        run_template(adam, dir = path$template_dir)

        # retrieve *.rda file in cache; copy to correct directory
        dataset_new <- load_rda(paste0(path$cache_dir, "/", adam, ".rda"))
        for (name in names(dataset_new)) {
            attr(dataset_new[[name]], "label") <- NULL
        }
        saveRDS(dataset_new, file = file.path("inst/verify/new", paste0(adam, ".rds")))

        
        dataset_old <- get_dataset_old(adam, path$adam_old_dir)
        # remove column attributes from old
        for (name in names(dataset_old)) {
            attr(dataset_old[[name]], "label") <- NULL
        }
        saveRDS(dataset_old, file = file.path("inst/verify/old", paste0(adam, ".rds")))

#        # temporary, for testing
         res = diffdf::diffdf(dataset_new, dataset_old, suppress_warnings = TRUE)
         if (diffdf::diffdf_has_issues(res)) {
            print(res)
            saveRDS(res, file=file.path("inst/verify", paste0(adam, ".diff")))
         }
 
   #     }, # end expr 
    #  error = function(e) {
    #       e
    #       message(paste0(adam , " problem: ", e$message))
    #      }
    #   )   # end tryCatch 

})  # end  purrr::map

    cli_inform("---- Done with templates\n")
} # end verify



#------------------------  helper functions

#' Display Results of running diffdf
#'
#' Reads all text files in a specified directory and returns their contents in a named list.
#' Each entry in the list contains the lines of one file as a character vector, and the
#' entry is named after the file.
#'
#' @param dir A character string specifying the path to the directory containing the text files.
#'
#' @return A named list. Each element is a character vector giving the lines of a file;
#'   the names of the list correspond to the individual filenames.
#'
#' @examples
#' \dontrun{
#' # Read contents of all files in "logs" directory
#' file_contents <- display_diff("logs")
#' print(names(file_contents)) # Shows the filenames
#' print(file_contents[[1]]) # Shows contents of the first file
#' }
#'
#' @param dir Directory with `diff` files, one for each ADaM
display_diff <- function(dir = NULL) {
  # Open connection to log file
  log_file <- file("verify_template_comparisons.txt", "a")

  files <- list.files(dir, full.names = TRUE)
  contents <- lapply(files, readLines)
  names(contents) <- basename(files)

  map2(
    names(contents), contents,
    function(name, content) {
      header <- paste("Differences found for ", str_replace_all(name, ".txt", ""),
                      " ", date(), "\n")

      # Display to console
      cli::cli_h1(header)
      cat(paste(content, collapse = "\n"))

      # Write to log file
      writeLines(header, log_file)
      writeLines(paste(content, collapse = "\n"), log_file)
      writeLines("\n", log_file) # Add separator between entries
    }
  )

  # Close the file connection
  writeLines("\n END of RUN -------", log_file)
  close(log_file)
  invisible(NULL)
}
#' @description: loads saved rda file `filename` and returns the dataset
load_rda <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"]) # returns dataset
}

#' @description saves the dataset to a specified file path
save_rda <- function(data, file_path, new_name) {
  # new_name must include  .rda
  if (missing(new_name)) {
    save(data, file = file_path, compress = "bzip2")
  }
}
#'
#' Loads and compares datasets using diffdf
#'
#' @desciption Uses `diffdf` to compare two datasets
#'  This package loads saved old rda files (generated from `get_dataset_old`) and uses
#'  diffdf to compare this data with a new dataset generated by `get_dataset_new`
#'  or a base for comparison.
#'
#' Usage
#'
#'   # Load two datasets using default comparison options
#'   get_dataset_new("path/to/adam.name.rda")
#'
#'   # Use custom comparison as file (e.g., with label files to be removed)
#'   compare(
#'     # The base for the comparison - a generic dataset object
#'     new_dataset = load_rda(path/to/base.object.rda),
#'
#'     # A file containing label information to remove from attributes of the datasets
#'     compare_file=get_R_data_path("_label_files/")
#' @param base name of the old ADaM dataset
#' @param compare name of the new ADaM dataset
#' @param keys set of keys to link the two datasets (optional)
#' @param file name of file to hold result of `diffdf`

compare <- function(base, compare, keys, file = NULL) {
  # DISCUSS ------------------------
  # Default: datasets are in tmp directories and not saved
  # To debug:  useful to have easy access to datasets in (1) global env
  e <- globalenv()
  e$old <- base
  e$new <- compare
  e$file <- file
  # ------------------------

  # remove column attributes
  for (name in names(base)) {
    attr(base[[name]], "label") <- NULL
  }

  }
#' Copies ADaM datasets from pharmaverseadam
#' @param adam_names character vector  Set of  ADaMs to download.
#' @param path Character string. Directory to save downloaded ADaMs.
download_adam_old <- function(adam_names, path = NULL) {

#  # NEW:
#  # ds here is singlular
#
#  f = function(ds) {
#      q = call("::", "pharmaverseadam", sym(ds))
#      save(q,
#          file = file.path(path, paste0(ds, ".rda"))) }
#
#  walk(.x = adam_names, .f = f)

#   LEGACY
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

}

#' Loads an ADaM dataset from a saved RDA file on disk.
#'
#' @param adam Character string. Name of the ADaM dataset to retrieve from
#'   disk (without file extension).
#' @param path Character string. Path to the directory containing the
#'   stored ADaM file.
#' @return The loaded dataset (as an R object).
#'
#' @examples
#' \dontrun{
#' adsl <- get_dataset_old("adsl", path = "data/adam")
#' }
get_dataset_old <- function(adam, path) {
  adam_old <- eval(load_rda(file.path(path, paste0(adam, ".rda"))))
}

#' Load a Saved ADaM Dataset

#' Loads an ADaM dataset from in  saved RDA file on disk.
#' @param adam Character string. Name of the ADaM dataset to retrieve
#'   from disk (without file extension).
#' @param path Character string. Path to the directory containing
#'   the stored ADaM file.
#'
#' @return The loaded dataset (as an R object).
get_dataset_new <- function(adam, path = NULL) {
  load_rda(paste0(path, "/", adam, ".rda"))
  adam_new <- get(adam)
}
#' Uses source() to run a template
#'
#' @param adam Character string.  Run the template associated with this ADaM name.
#' @param dir Character string.  Directory where templates are stored.
#' @return results of running template
#' @description  runs the template for the specified ADaM
run_template <- function(adam, dir = NULL) {
  source(paste0(dir, "/ad_", adam, ".R")) # nolint
}

# Since this is script, and not package R function it must be loaded/sourced separately.
# The next line runs the entire process after this file is sourced.
# verify_templates()   # nolint


#' create_directories
#'
#' @return list of paths to directories
create_directories <- function() {
  cli_inform("Creating temporary directories")
  x <- tempdir()
  if (!dir.exists(file.path(x, "old") )) dir.create(file.path(x, "old"), showWarnings = TRUE)
  if (!dir.exists(file.path(x, "new"))) dir.create(file.path(x, "new"), showWarnings = TRUE)
  if (!dir.exists(file.path(x, "diff"))) dir.create(file.path(x, "diff"), showWarnings = TRUE)

  list(
    template_dir = "inst/templates",
    cache_dir = tools::R_user_dir("admiral_templates_data", which = "cache"),
    adam_new_dir = file.path(x, "new"),
    adam_old_dir = file.path(x, "old"),
    diff = file.path(x, "diff")
  )
}
    
