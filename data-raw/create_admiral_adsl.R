#  Create dataset:   data/admiral_adsl.rda
#  This script, create_admiral_adsl.R, creates dataset data/admiral_adsl.rda.
#

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
