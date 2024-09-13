#  Create dataset:   data/admiral_adlb.rda


# The ultimate goal of this script (`create_admiral_adlb.R`) is to create dataset data/admiral_adlb.rda.
# This script controls the multi-step process.


# First, use template to create the R script (in data-raw/admiral_adlb.R).
# Next, source this script and create the data (~/.cache/R/admiral_template_data/admiral_adlb.rda)
# Finally, shorten this data (now ~ 1.2 MB) by selecting only certain USUBJID


# First,  create the R script (from a template)
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

# Load the data into .GlobalEnv
load("~/.cache/R/admiral_templates_data/adlb.rda")

# nrow(adlb)   83,652

# limit rows, by selecting only these USUBJID
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

# check
# USUBJID %in% result$USUBJID

# Finally, save reduced ds
use_data(admiral_adlb, overwrite = TRUE)

##
##  TEST - is dataset identical to .... backup of unaltered dataset
##
e1 = new.env()
e2 = new.env()

load("data/admiral_adlb.rda", e1)

# CHANGE to your location of original dataset
load("data-backup/admiral_adlb.rda", e2)

identical(e1$admiral_adsl, e2$admiral_adsl)

# cleanup
rm(e1)
rm(e2)
