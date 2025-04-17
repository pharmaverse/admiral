# This script:  data-raw/admiral_verify_templates.R

#
# Generates ADaM from templates and compares to previously generated ADaM file (found in pharmaverseadam).
# pharmaverseadam is the SOURCE.
#
# Much code taken from pharamavreseadam::create_adams_data.R
# (https://github.com/pharmaverse/pharmaverseadam/blob/main/data-raw/create_adams_data.R)
#
# USAGE:   verify_templates()
#
# TODO:
# - add `attr` to generated ADaM
# - re-write `load_rda()`:  now serves two purposes and was written in simpler times
# - where code overlaps, make pharamversadam script and this one the same.

# Assumptions/Questions:
# - ignore *.rda files in admiral/data (per Ben)
# - written as standalone R script, not as package R function
# - compares full ADaM - all rows (ie no reduction in number of rows in each dataset)
# - use cli:: for messages/errors?

#' (if were to add to `admiral` package)
#' @param pkg  package (ex:  "admiral )
#' @param adams_names ADaM or CDISC name, without prefix or suffix  (ex:  adlb)
#' @param adams_new  character vector of ADaM after template is run
#' @param adams_old  character vector of original ADaM done at earlier date, saved in github

#' @param template_dir path  Directory where templates where package is installed (not sourced)
#' @param cache_dir (cache)     ~/.config/R/admiral_templates_data   (varies by OS/configuration)

#' Directories
#' - admiral/data     not used
#' - admiral/inst/templates/  templates to create ADaM datasets
#' - cache_dir  cache, where newly generated ADaM files kept and diffdf reports
#'      (os dependent)
#' - adam_old_dir  - dir where prior ADaMs are downloaded and stored

load_rda <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"]) # returns dataset
}
save_rda <- function(data, file_path, new_name) {
  # new_name must include  .rda
  if (missing(new_name)) {
    save(data, file = file_path, compress = "bzip2")
  }
}
compare <- function(base, compare, keys, file=NULL) {
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

clean_adam_old_dir = function(){
  if (dir.exists(tempdir())) {
    unlink(adam_old_dir, recursive = TRUE)
    message(adam_old_dir,  "directory deleted: ", tempdir())
  } else {
    message(adam_old_dir, "directory does not exist: ", tempdir())
  }
}

download_adam_old = function (adam_names){
  # tempdir = directory to store
  # adam_names = chr[] of old adams to download
  # tempdir() is fixed for session
  # DO NOT LOAD to memory these files till needed.

  lapply(adam_names, function(adam){

    githubURL = paste0(
      "https://github.com/pharmaverse/pharmaverseadam/raw/refs/heads/main/data/",
      adam, ".rda?raw=true")

    cat(adam, "\n")
    download.file(url = githubURL,
                  quiet = TRUE,
                  destfile = paste0(path$adam_old_dir, "/", adam, ".rda"),
                  mode = "wb")
  })
}
  # does the real work
  run_template <- function(adam) {
    source(paste0(path$template_dir, "/ad_", adam, ".R"))
  }

get_dataset_old = function(adam) {
  # TODO: replace with LATEST files from github, not package.
  #  adam_old <- do.call(`::`, args = list("pharmaverseadam", adam))
  adam_old <- load_rda(paste0(path$adam_old_dir, "/", adam, ".rda"))
}

get_dataset_new = function(adam){
    load_rda(paste0(path$cache_dir, "/", adam, ".rda"))
    adam_new <- get(adam)
}

verify_templates <- function(pkg = "admiral", ignore_templates_pkg = NULL) {

  #   SETUP ----
  # TODO: remove prior ADaM downloads
  clean_cache()

  #   TODO ----
  # DISCUSS:  Diffdf creates *.txt files, DIRECTORY to put them?
  # TODO:   fct to remove remove old diffdf *.txt files

  # testing ----
  pkg = "admiral"

  # DISCUSS:  more pkg checking?
  if (pkg != "admiral") error("Curently, only admiral package is accepted.")

  library(pkg, character.only = TRUE)
  library(teal.data)
  sprintf("generating ADaMs for  %s package\n", pkg)

  # temporary directories
  x = tempdir()
  dir.create(paste0(x, "/old"))
  dir.create(paste0(x, "/new"))
  dir.create(paste0(x, "/diff"))

  # list of important paths
  path <<- list(
    template_dir = file.path(system.file(package = pkg), "templates"),
    cache_dir = tools::R_user_dir("admiral_templates_data", which = "cache"),
    adam_old_dir =  paste0(x, "/old"),
    adam_new_dir =  paste0(x, "/new"),
    diff = paste0(x, "/diff")
  )

  # gather all templates for this pkg (12 found) ----
  templates <- list.files(path$template_dir, pattern = "ad_")
  templates


  #adam_names is unnamed chr[]
  ## from templates  generate these ADaMs, 12 ----
  adam_names <- vapply(templates, function(x) gsub("ad_|\\.R", "", x), USE.NAMES = FALSE, character(length = 1))
  adam_names

  # check
  if (length(templates) != length(adam_names))  stop("Number of templates and adam_names differ")

  # templates is a named chr[]
  names(templates) <- adam_names
  templates


  # per Ben, ignore "adlbhy"
  adam_names = adam_names[adam_names != "adlbhy"]
  adam_names

  # download, save prior ADaMs from pharmaverseadam
  download_adam_old(adam_names)

  # keys for diffdf
  keys <- teal.data::default_cdisc_join_keys

  # keys is named list of keys, each element (for each adam_name) and is chr[] of keys
  keys =sapply(toupper(adam_names), function(e) unname(keys[[e]][[e]]), USE.NAMES=T)
  names(keys) <- tolower(names(keys))
  keys


  # finally, construct a named list object, each element (adam_names) holds a template and keys
  obj=lapply(adam_names, function(e) list("template" = templates[[e]],
    "keys" = keys[[e]]))
  names(obj) = adam_names

  # done !
  obj
  names(obj)


  # testing check ----
  if (F) {
  # need info for specific adam_name?
  obj[["adae"]]$template
  obj[["adae"]]$keys
  }

  # DISCUSS:   NOW, user must manually enter templates to ignore ----
  # CHANGE TO: PROPOSED:   user must manually remove ADaMs to ignore

  # CHANGE:   remove this -----------
  if (F) {
  # FOR TESTING, (omit ad_adlb.R)
  # ignore_templates_pkg = templates[!(templates %in% c("ad_adsl.R", "ad_adlb.R"))]
  # For real,
  ignore_templates_pkg <- templates[templates == c("ad_adlbhy.R")]
  }
  # CHANGE:   remove above -----------


  # CHANGE:  now, which templates to run?   ENTER in adam_names
  # (for TESTING, skip "adlb" is lengthy ----)

  # for TESTING, otherwise will be ALL adam_names
  adam_names = c("adsl", "adae")

  sprintf("---- Run templates\n")
  compare_list <- purrr::map(adam_names, .progress = TRUE, function(adam){
    cat("---- Template running for ", adam, "\n")
    run_template(adam)
    dataset_new <- get_dataset_new(adam) # this function would retrieve that dataset from cache
    dataset_old <- get_dataset_old(adam)
    compare(base=dataset_old, compare = dataset_new, keys=obj[[adam]]$keys,
            file=paste0(path$diff,"/", adam,".txt"))
    # TODO: remove old/new ADaMs from environment
  })

  print("DONE")
}

