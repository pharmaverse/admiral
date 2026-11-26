#  Create dataset:   data/admiral_adsl.rda
#  This script, create_admiral_adsl.R, creates dataset data/admiral_adsl.rda.
#

# Create dataset data/admiral_adsl.rda

# Run template script to create adsl
source("inst/templates/ad_adsl.R", echo = TRUE) # nolint

admiral_adsl <- adsl # use correct name

# Store the metadata, so that the dataset prints and summarizes as an ADSL
admiral_adsl <- admiral_df(
  admiral_adsl,
  dataset_name = "ADSL",
  dataset_class = "ADSL",
  dataset_keys = exprs(STUDYID, USUBJID)
)

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
