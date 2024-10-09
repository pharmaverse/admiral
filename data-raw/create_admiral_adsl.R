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


# Create dataset data/admiral_adsl.rda

# Run template script to create adsl
source("inst/templates/ad_adsl.R", echo = TRUE) # nolint

admiral_adsl <- adsl # use correct name

# Get previous dataset for comparison
adsl_old <- admiral::admiral_adsl

# Finally, save reduced dataset
usethis::use_data(admiral_adsl, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = adsl_old,
  compare = admiral_adsl,
  keys = c("STUDYID", "USUBJID")
)
