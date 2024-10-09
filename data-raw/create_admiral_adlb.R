#  Create dataset:   data/admiral_adlb.rda
#  This script:  create_admiral_adlb.R creates dataset data/admiral_adlb.rda.
#

# Preliminary
library(diffdf) # nolint

# To clarify directories (can be removed)
# nolint start: object_name_linter
CACHE_DIR <- tools::R_user_dir("admiral_templates_data", which = "cache")
DATA_DIR <- "data-dir"
DATA_RAW <- "data-raw"
TEMPLATE_DIR <- "inst/templates/"

# clean CACHE_DIR
THE_FILE <- paste0(CACHE_DIR, "/adlb.rda")
THE_FILE <- paste0(CACHE_DIR, "/adsl.rda")
if (file.exists(THE_FILE)) file.remove(THE_FILE)
# nolint end


# Create dataset data/admiral_adlb.rda

# Run template script to create adlb
source("inst/templates/ad_adlb.R", echo = TRUE) # nolint

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

admiral_adlb <- filter(adlb, USUBJID %in% usubjids)


# Get previous dataset for comparison
adlb_old <- admiral::admiral_adlb

# Finally, save reduced dataset
usethis::use_data(admiral_adlb, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = adlb_old,
  compare = admiral_adlb,
  keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT")
)
