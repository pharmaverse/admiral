# This script:  data-raw/admiral_verify_templates.R
#
# Generates ADaM from templates and compares to previously generated ADaM in pharmaverseadam.
# pharmaverseadam is the SOURCE.
#
# Much code taken from pharamavreseadam::create_adams_data.R
# (https://github.com/pharmaverse/pharmaverseadam/blob/main/data-raw/create_adams_data.R)
#

# USAGE:   verify_templates()
#
# TODO:
# - add `attr` to generated ADaM
# - where code overlaps, make pharamversadam script and this one the same.
# - vectorize, using lapply

# Assumptions/Questions:
# - ignore *.rda files in admiral/data (per Ben)
# - written as standalone R script, not as package R function
# - compares full ADaM - all rows (ie no reduction in number of rows in each dataset)
# - change ?   in template ad_adsl.R, change variable `dir`?
# - now using `.rda` datafiles, switch to `.rds`? (easier to program; pilot5 uses)
# - use cli:: for messages/errors?

#' (if added to `admiral` package)
#' @param pkg  package (ex:  "admiral )
#' @param name ADaM or CDISC name, without prefix or suffix  (ex:  adlb)
#' @param adam_new  ADaM after template is run
#' @param adam_old =  original ADaM done at earlier date, saved in pharamverseadam

#' @param tp   template R file  (ex: ad_adlb.R)
#' @param templates_path path  Directory where templates where package is installed (not sourced)
#' @param dataset_dir (cache)     ~/.config/R/admiral_templates_data   (varies by OS/configuration)

#' Directories
#' - admiral/data     not used
#' - admiral/inst/templates/  templates to create ADaM datasets
#' - dataset_dir  cache, where newly generated ADaM files kept and diffdf reports
#'      (os dependent)

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
compare <- function(base, compare, keys, file) {
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
  dataset_dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
  if (dir.exists(dataset_dir)) {
    unlink(dataset_dir, recursive = TRUE)
    message("Cache directory deleted: ", dataset_dir)
  } else {
    message("Cache directory does not exist: ", dataset_dir)
  }
}

#  IGNORE:  not in use
reduce <- function() {
  ## DISCUSSION:   adlb.rda ADaM is large and must be reduced in size (for testing)
  ## What about other ADaMs ?  use vector
  if (tp == "ad_adlb.R") {
    # Reduce size of adlb dataset for testing
    # Limit rows by selecting only these USUBJIDs
    usubjids <-
      c(
        "01-701-1015",
        "01-701-1023",
        "01-701-1028",
        "01-701-1033",
        "01-701-1034",
        "01-701-1047",
        "01-701-1097",
        "01-705-1186",
        "01-705-1292",
        "01-705-1310",
        "01-708-1286"
      )
    # small adlb
    admiral_adlb <- dplyr::filter(adlb, USUBJID %in% usubjids)
    adam_new <- admiral_adlb
  } else {
    admiral_adsl <- adsl
    adam_new <- admiral_adsl
  } ## end if
} ## end reduce

verify_templates <- function(pkg = "admiral", ignore_templates_pkg = NULL) {
  clean_cache()
  library(pharmaverseadam)

  # TODO:  more checking
  if (pkg != "admiral") error("Curently, only admiral package is accepted.")

  # SOURCE:  Use pharmaverseadam as source for ADaM *.rda files (22 ADaM files)
  # ignore adlbhy.rda (per Ben)
  sprintf("generating ADaMs for  %s package", pkg)
  library(pkg, character.only = TRUE)

  source_adams <- data(package = "pharmaverseadam")

  # 23 ADaMs found
  source_adams <- source_adams$results[, "Item"]

  # per Ben, ignore "adlbhy"
  # source_adams = source_adams[source_adams != "adlbhy"]
  source_adams

  # gather keys (for diffdf)
  # DISCUSS use of keys in CDISC
  keys_adsl <- teal.data::default_cdisc_join_keys[c("ADSL")]$ADSL
  keys_adsl

  # gather all templates for this pkg (12 found)
  template_path <- file.path(system.file(package = pkg), "templates")
  templates <- list.files(template_path, pattern = "ad_")
  templates

  ## BUT, |templates|  <  |ADaMs in pharmaverseadam|   !
  ## Match/Reduce number of source_adams !

  ## templates can generate these ADaMs, 12
  names <- vapply(templates, function(x) gsub("ad_|\\.R", "", x), USE.NAMES = FALSE, character(length = 1))

  # KEEP only ADaMs with template (now 12)
  source_adams <- source_adams[source_adams %in% names]

  #
  # finally name (ADaM or CDISC name)  the templates vector
  names(templates) <- names


  # list of important paths
  path <- list(
    "templates_path" = file.path(system.file(package = pkg), "templates"),
    dataset_dir = tools::R_user_dir("admiral_templates_data", which = "cache")
  )

  # For testing, (omit ad_adlb.R)
  # ignore_templates_pkg = templates[!(templates %in% c("ad_adsl.R", "ad_adlb.R"))]

  # For real,
  ignore_templates_pkg <- templates[templates == c("ad_adlbhy.R")]

  #
  # now, generate ADaM datasets and put in cache dir
  #

  # begin tp  ----
  for (tp in templates) {
    # Each tp creates single ADaM package
    # if (tp != "ad_adsl.R") next
    # get CDISC name
    name <- gsub("ad_|\\.R", "", tp)
    adam_old <- do.call(`::`, args = list("pharmaverseadam", name))
    # run template, which caches new adam_new
    if (tp %in% ignore_templates_pkg) {
      next
    } else {
      cat("Running template: ", tp, "\n")
      source(file.path(template_path, tp), echo = TRUE) # nolint
    }
    # retrieve new adam from cache dir (puts into globalenv())
    load_rda(paste0(path$dataset_dir, "/", name, ".rda"))
    adam_new <- get(name)

    ## TODO function use vectors
    sprintf("comparing ... diffdf")
    cat("comparing ...", tp, "\n")
    #   compare(base=adam_old,
    #            compare=adam_new,
    #           keys= ifelse(tp == "adlb.R",
    #                        c("USUBJID", "PARAMCD", "AVISIT", "ADT"),
    #                        c("STUDYID","USUBJID")),
    #            file = paste0("data/diff_", name, ".txt")
    diffdf::diffdf(
      base = adam_old,
      compare = adam_new,
      keys = NULL,
      file = paste0("data/diff_", name, ".txt")
    )
  } ## end tp

  print("DONE")
} ## end main
