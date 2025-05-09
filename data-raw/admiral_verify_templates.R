# This script:  data-raw/admiral_verify_templates.R

# FIX
# Generates ADaM from templates and compares to previously generated ADaM file
# (These are found in https:://github.com main branch pharmaverseadam).
#
#
# Much code taken from pharamavreseadam::create_adams_data.R
# (https://github.com/pharmaverse/pharmaverseadam/blob/main/data-raw/create_adams_data.R)
#
# USAGE:   verify_templates()
#
# TODO:
# - add `attr` to generated ADaM
# - where code overlaps, make pharamversadam script and this one the same.

# Assumptions/Questions:
# - ignore *.rda files in admiral/data (per Ben)
# - compares full ADaM - all rows (ie no reduction in number of rows in each dataset)
# - use cli:: for messages/errors - YES
# - use has loaded/attached `admiral` (library(admiral)). The script is NOT being run independently.

#' (if were to add to `admiral` package)
#' @param pkg  package (ex:  "admiral )
#' @param adams_names ADaM or CDISC name, without prefix or suffix  (ex:  adlb)
#' @param adams_new  character vector of ADaM after template is run
#' @param adams_old  character vector of original ADaM done at earlier date, saved in github

## DISCUSS
#' @param template_dir_path  installed directory (admiral/inst/templates/)
#' @param cache_dir (cache)     ~/.config/R/admiral_templates_data   (varies by OS/configuration)
#' @param adam_old_dir  temporary directory where old ADaMs (downloaded from github.com) are stored
#' @param adam_new_dir  temporary directory where new ADaMs are stored.

verify_templates <- function(pkg = "admiral", ignore_templates_pkg = NULL) {
  # SETUP ----
  # TODO: remove prior ADaM downloads
  # ASSUME:  (1) user is running script for 1st time and no temporary directories exist yet, OR
  #          (2) user is running script a 2nd time and must clear old

  clean_cache() # clear all..

  # TODO:   fct to remove remove old diffdf *.txt files

  pkg <- "admiral"
  if (pkg != "admiral") error("Curently, only admiral package is accepted.")

  library(pkg, character.only = TRUE)
  library(teal.data)
  library(purrr)
  sprintf("generating ADaMs for  %s package\n", pkg)

  # temporary directories
  x <- tempdir()
  dir.create(paste0(x, "/old"))
  dir.create(paste0(x, "/diff"))

  # TODO: choose 1:
  # cache_dir and adam_new_dir are the SAME
  # cache_dir is where templates deposit new ADaM
  # adam_new_dir is used when reading from disk
  # list of important paths
  ## path <<- list(
  path <- list(
    template_dir = file.path(system.file(package = pkg), "templates"),
    cache_dir = tools::R_user_dir("admiral_templates_data", which = "cache"),
    adam_new_dir = tools::R_user_dir("admiral_templates_data", which = "cache"),
    adam_old_dir = paste0(x, "/old"),
    diff = paste0(x, "/diff")
  )

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

  # per Ben, ignore "adlbhy"
  adam_names <- adam_names[adam_names != "adlbhy"]

  # download, save prior ADaMs from pharmaverseadam

  download_adam_old(adam_names, path = path$adam_old_dir)

  # keys for diffdf
  keys <- teal.data::default_cdisc_join_keys

  # keys is named list of keys, each element (for each adam_name) and is chr[] of keys
  keys <- sapply(toupper(adam_names), function(e) unname(keys[[e]][[e]]), USE.NAMES = TRUE)
  names(keys) <- tolower(names(keys))
  keys


  # finally, construct a named list object: each element (adam_names) holds a template and keys
  obj <- lapply(adam_names, function(e) {
    list(
      "template" = templates[[e]],
      "keys" = keys[[e]]
    )
  })
  names(obj) <- adam_names

  # done !
  obj
  names(obj)

  # for TESTING, otherwise will be ALL adam_names
  adam_names <- c("adsl", "adae")

  sprintf("---- Run templates\n")
  compare_list <- purrr::map(adam_names, .progress = TRUE, function(adam) {
    cat("---- Template running for ", adam, "\n")
    run_template(adam, dir = path$template_dir)
    dataset_new <- get_dataset_new(adam, path$cache_dir)
    dataset_old <- get_dataset_old(adam, path$adam_old_dir)
    compare(
      base = dataset_old,
      compare = dataset_new,
      keys = obj[[adam]]$keys,
      file = paste0(path$diff, "/", adam, ".txt")
    )
  })

  # TODO:  cleanup?

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
#' @export
display_diff <- function(dir = NULL) {
  files <- list.files(dir, full.names = TRUE)
  contents <- lapply(files, readLines)
  names(contents) <- basename(files)
  contents
}

load_rda <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"]) # returns dataset
}
save_rda <- function(data, file_path, new_name) {
  # new_name must include  .rda
  if (missing(new_name)) {
    save(data, file = file_path, compress = "bzip2")
  }
}
compare <- function(base, compare, keys, file = NULL) {
  tryCatch(
    {
      diffdf::diffdf(
        base = base,
        compare = compare,
        keys = keys,
        file = file
      )
    },
    error = function(e) {
      message("Error in diffdf: ", e$message)
    }
  ) ## end tryCatch
} ## end compare

# start fresh
clean_cache <- function() {
  cache_dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
  if (dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    message("Cache directory deleted: ", cache_dir)
  } else {
    message("Cache directory does not exist: ", cache_dir)
  }
}

clean_adam_old_dir <- function(dir = NULL) {
  if (dir.exists(tempdir())) {
    unlink(adam_old_dir, recursive = TRUE)
    message(adam_old_dir, "directory deleted: ", tempdir())
  } else {
    message(adam_old_dir, "directory does not exist: ", tempdir())
  }
}

download_adam_old <- function(adam_names, path = NULL) {
  lapply(adam_names, function(adam) {
    githubURL <- paste0(
      "https://github.com/pharmaverse/pharmaverseadam/raw/refs/heads/main/data/",
      adam, ".rda?raw=true"
    )
    cat("Downloading adam from github pharamaversedam", adam, "\n")
    download.file(
      url = githubURL,
      quiet = TRUE,
      destfile = paste0(path, "/", adam, ".rda"),
      mode = "wb"
    )
  })
}

get_dataset_old <- function(adam, path = NULL) {
  adam_old <- load_rda(paste0(path, "/", adam, ".rda"))
}

get_dataset_new <- function(adam, path = NULL) {
  load_rda(paste0(path, "/", adam, ".rda"))
  adam_new <- get(adam)
}

# does the real work
run_template <- function(adam, dir = NULL) {
  source(paste0(dir, "/ad_", adam, ".R"))
}
