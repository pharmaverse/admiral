#  Create dataset:   data/admiral_adlb.rda
#  This script:  create_admiral_adlb.R creates dataset data/admiral_adlb.rda.
#

# Preliminary
library(diffdf) # nolint

# To clarify directories (can be removed)
# nolint start: object_name_linter
CACHE_DIR <- "~/.cache/R/admiral_templates_data/"
DATA_DIR <- "data-dir"
DATA_RAW <- "data-raw"
TEMPLATE_DIR <- "inst/templates/"

# clean CACHE_DIR
THE_FILE <- paste0(CACHE_DIR, "/adlb.rda")
THE_FILE <- paste0(CACHE_DIR, "/adsl.rda")
if (file.exists(THE_FILE)) file.remove(THE_FILE)
# nolint end

#
# STEPS
#
# First, use template to create the R script (in data-raw/admiral_adlb.R).
# Next, source this script and create the data (~/.cache/R/admiral_template_data/admiral_adlb.rda)
# Finally, shorten this data (now ~ 1.2 MB) by selecting only certain USUBJID

# orignal method  - OMIT
if (FALSE) {
  # ### original mehtod (method 1)
  # # First,  create the R script (from a template)
  adam_name <- "adlb"
  save_path <- paste0("./data-raw/admiral_", adam_name, ".R")
  use_ad_template(
    adam_name = adam_name,
    save_path = save_path,
    open = FALSE,
    overwrite = TRUE
  )
  # Second, source the script and save data in .cache
  source("data-raw/admiral_adlb.R") # nolint
  load("~/.cache/R/admiral_templates_data/adlb.rda")
}

#
# Instead, USE template, as recommened by Buzz
#
source(paste0(TEMPLATE_DIR, "/ad_adlb.R")) # nolint
load(paste0(CACHE_DIR, "adlb.rda"))

#
# limit rows, by selecting only these USUBJID
#
#' 01-701-1015, 01-701-1023, 01-701-1028, 01-701-1033,
#' 01-701-1034, 01-701-1047, 01-701-1097, 01-705-1186,
#' 01-705-1292, 01-705-1310, 01-708-1286

usubjid <-
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

#  prepare for inner join
user <- tibble(
  USUBJID = usubjid
)
result <- inner_join(adlb, user)
admiral_adlb <- result
admiral_adlb

#
# Finally, save reduced ds
#
usethis::use_data(admiral_adlb, overwrite = TRUE)
