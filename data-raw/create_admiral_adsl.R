#  Create dataset:   data/admiral_adsl.rda
#  This script, create_admiral_adsl.R, creates dataset data/admiral_adsl.rda.
#
# Preliminary
library(diffdf) # nolint
#
# nolint start: object_name_linter.
# To clarify directories (can be removed)
CACHE_DIR <- "~/.config/cache/R/admiral_templates_data/"
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
#
source(paste0(TEMPLATE_DIR, "/ad_adsl.R")) # nolint
load(paste0(CACHE_DIR, "adsl.rda"))
admiral_adsl <- adsl

#
# Finally, save reduced ds
#
usethis::use_data(admiral_adsl, overwrite = TRUE)
