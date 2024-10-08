#  Create dataset:   data/admiral_adsl.rda
#  This script, create_admiral_adsl.R, creates dataset data/admiral_adsl.rda.
#
# Preliminary
library(diffdf) # nolint
#
# nolint start: object_name_linter.
# To clarify directories (can be removed)
CACHE_DIR <- tools::R_user_dir("admiral_templates_data", which = "cache")
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
# Source template (ad_adsl.R) to create the data (adsl.rda) and save in CACHE_DIR
#
source(paste0(TEMPLATE_DIR, "/ad_adsl.R")) # nolint
load(paste0(CACHE_DIR, "/adsl.rda"))
admiral_adsl <- adsl # use correct name

#
# Finally, save reduced ds
#
usethis::use_data(admiral_adsl, overwrite = TRUE)
