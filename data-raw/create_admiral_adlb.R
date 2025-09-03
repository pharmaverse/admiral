#  Create dataset:   data/admiral_adlb.rda
#  This script:  data-raw/create_admiral_adlb

# Create dataset data/admiral_adlb.rda

# First, retrieve adlb dataset as it is now
adlb_old <- admiral::admiral_adlb

# Run template script to create new adlb.rda
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

admiral_adlb <- dplyr::filter(adlb, USUBJID %in% usubjids)

# Finally, save reduced dataset
usethis::use_data(admiral_adlb, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = adlb_old,
  compare = admiral_adlb,
  keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT")
)
