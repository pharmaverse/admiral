#  Create dataset:   data/admiral_adlb.rda
#  This script:  data-raw/create_admiral_adlb



# adlb_old  adlb as it in admiral/data before do anything
# adlb      after template is run
# admiral_adlb after template is run and reduced
# cache     ~/.config/R/admiral_templates_data   (varies by OS/configuration)
# admiral/data     dir where admiral_adlb.rda will be saved
# admiral/inst/templates/ad_adlb.R  template to create adlb

#
# Create dataset data/admiral_adlb.rda

# First, retrieve adlb dataset as it is now
adlb_old <- admiral::admiral_adlb

# Next, run template script to create new adlb.rda/ save in cache/ put adlb into globalenv()
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


# Finally, save reduced dataset in data/
usethis::use_data(admiral_adlb, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = adlb_old,
  compare = admiral_adlb,
  keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT")
)

