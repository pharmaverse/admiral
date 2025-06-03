# This script:  data-raw/admiral_verify_templates.R

# TODO:
# - where code overlaps, make pharamversadam script and this one the same.
# - TODO: replace messages with cli::cli_*
# - TODO: lots of directories and named lists.  Best way to keep neat and organized?

# Assumptions/Questions:
# - ignore *.rda files in admiral/data (per Ben)
# - compares full ADaM - ie all rows
# - for developer use and developer has run load_all()
# - use cli:: for messages/errors - YES
# - most datasets are  stored as.rda files


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

## DISCUSS, several named vectors, try tibble?
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
#' @param ds  character vector of ADaM names.  Corresponding to templates to run.  \
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
verify_templates <- function(pkg = "admiral", ds = c("adae")) {
  # TODO: delete all remove prior ADaM downloads
  # ASSUME:  (1) user is running script for 1st time and no temporary directories exist, OR
  #          (2) user is running script a 2nd time, in same session, and must remove directories

  clean_cache() # clear all..

  pkg <- "admiral"
  if (pkg != "admiral") cli_abort("Currently only `admiral` package is supported")

  # nolint start
  library(pkg, character.only = TRUE)
  library(teal.data)
  library(purrr)
  library(cli)
  # nolint end
  cli_alert("Generating ADaMs for { pkg} package.")

  # TODO:  do not show Warnings
  # (DISCUSS) CHOOSE one:  create temporary or permanent directories
  cli_inform("Creating temporary directories")
  x <- tempdir()
  #x  <- "~/data"
  dir.create(file.path(x, "old"), showWarnings = TRUE)
  dir.create(file.path(x, "new"), showWarnings = TRUE)
  dir.create(file.path(x, "diff"), showWarnings = TRUE)

  # TODO:
  path <- list(
    template_dir = "inst/templates",
    cache_dir = tools::R_user_dir("admiral_templates_data", which = "cache"),
    adam_new_dir = file.path(x, "new"),
    adam_old_dir = file.path(x, "old"),
    diff = file.path(x, "diff")
  )

  # TODO:
  # if dir exists then empty it ; if not exist then create it.
  # lapply(path, function(x)  if( x %in% c("template_dir", "cache_dir")!exists(x)) dir.create(x))

  # gather all templates for this pkg (12 found) ----
  templates <- list.files(path$template_dir, pattern = "ad_")

  # from templates  generate vector of adam_names
  adam_names <- vapply(templates, function(x) gsub("ad_|\\.R", "", x),
    USE.NAMES = FALSE, character(length = 1)
  )

  # check
  if (length(templates) != length(adam_names)) stop("Number of templates and adam_names differ")

  # templates is a named chr[]
  names(templates) <- adam_names

  # User specified templates in `ds`. Limit to only these templates.
  adam_names <- adam_names[adam_names %in% ds]

  # per Ben, ignore "adlbhy"
  adam_names <- adam_names[adam_names != "adlbhy"]

  # download, save prior ADaMs from pharmaverseadam as .rda files
  cli_inform("---- Begin downloading from github pharmaverseadam")
  download_adam_old(adam_names, path = path$adam_old_dir)

  cli_inform("---- Run templates\n")
  compare_list <- purrr::map(adam_names, .progress = TRUE, function(adam) {
    cli_inform("Template running for {adam}")
    run_template(adam, dir = path$template_dir)
    # retrieve *.rda file in cache; copy to correct directory
    dataset_new = load_rda(paste0(path$cache_dir, "/", adam, ".rda"))
    file.copy(file.path(path$cache_dir, paste0(adam, ".rda")),
              file.path(path$adam_new_dir, paste0(adam, ".rda")))
    dataset_old <- get_dataset_old(adam, path$adam_old_dir)

    # ------------------------  DISCUSS
    # CHOOSE ONE:   continue to compare,  
    # OR,  insert Eli's quarto code to compare AND display nicely
    # ------------------------  
    compare(
      base = dataset_old,
      compare = dataset_new,
      keys = teal.data::default_cdisc_join_keys[[adam]],
      file = paste0(path$diff, "/", adam, ".txt")
    )
  })
  # run AFTER verify_templates() completes, otherwise will not print
  display_diff(dir = path$diff)
}

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
display_diff <- function(dir = NULL) {
  files <- list.files(dir, full.names = TRUE)
  contents <- lapply(files, readLines)
  names(contents) <- basename(files)
  map2(
    names(contents), contents,
    function(name, content) {
      cli::cli_h1(paste(
        "Differences found for", str_replace_all(name, ".txt", ""),
        " ", today(), "\n"
      ))
      cat(paste(content, collapse = "\n"))
    }
  )
  invisible(NULL)
}

#' @description: loads saved file `filename` and returns the dataset
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
#'  This package loads saved RDA files (generated from `get_dataset_old`) and uses
#'  diffdf to compare the loaded data with a new dataset generated by `get_dataset_new`
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
compare <- function(base, compare, keys, file = NULL) {
  # DISCUSS
  #--------debugging--------- remove OR save as new save_debug( , debug=FALSE) function
  # useful to have easy access to these in global env
  e <- globalenv()
  e$old <- base
  e$new <- compare
  e$file <- file
  # nolint start
  # saveRDS(e$old,file= paste0("old", ".RDS"))    # temporary
  # saveRDS(e$new,file= paste0("new", ".RDS"))
  # ------------------------  to be removed ?
  # nolint end

  # remove column attributes
  for (name in names(base)) {
    attr(base[[name]], "label") <- NULL
  }
  for (name in names(compare)) {
    attr(compare[[name]], "label") <- NULL
  }
  tryCatch(
    {
      e$res <- diffdf::diffdf(
        base = base,
        compare = compare,
        keys = keys,
        file = file,
        suppress_warnings = TRUE # for now
      )
    },
    error = function(e) message("Error in diffdf: ", e$message)
  ) ## end tryCatch
} ## end compare

#' @description:  removes the cache directory
clean_cache <- function() {
  cache_dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
  if (dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    message("Cache directory deleted: ", cache_dir)
  } else {
    message("Cache directory does not exist: ", cache_dir)
  }
}

#' @description:  removes the old adam directory
clean_adam_old_dir <- function(dir = NULL) {
  if (dir.exists(tempdir())) {
    unlink(adam_old_dir, recursive = TRUE)
    message(adam_old_dir, "directory deleted: ", tempdir())
  } else {
    message(adam_old_dir, "directory does not exist: ", tempdir())
  }
}

#' Downloads ADaM datasets from pharmaverseadam
#' @param adam_names character vector  Set of  ADaMs to download.
#' @param path Character string. Directory to save downloaded ADaMs.
download_adam_old <- function(adam_names, path = NULL) {
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

## DISCUSS:
#  Combine next 2 functions into one.   get_dataset(adam, path)
#  The path tells us if old or new.

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
  adam_old <- load_rda(file.path(path, paste0(adam, ".rda")))
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
#' @param adam Character string.  Run the template assicated with this ADaM name.
#' @param dir Character string.  Directory where templates are stored.
#' @return
#' @description  runs the template for the specified ADaM
run_template <- function(adam, dir = NULL) {
  source(paste0(dir, "/ad_", adam, ".R")) # nolint
}

# Since this is script, and not package R function it must be loaded/sourced separately.
# The next line runs the entire process after this file is sourced.
#verify_templates()   # nolint
