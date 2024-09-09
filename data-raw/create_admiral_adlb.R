#  Create dataset:   data/admiral_adlb.rda

# This is a MULTI-step process.

# First, using template to create the R script (in data-raw/admiral_adlb.R) which will generate the data.
# Next, source this script and create the data (~/.cache/R/admiral_template_data/admiral_adlb.rda)
# Finally, shorten this data (now ~ 1.2 MB) by selecting only certain USERJID


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
source("data-raw/admiral_adlb.R")  # nolint

# Load the data into .GlobalEnv
load("~/.cache/R/admiral_templates_data/adlb.rda" )

#nrow(adlb)   83,652

# limit rows, by selecting only these USUBJID
#' 01-701-1015, 01-701-1023, 01-701-1028, 01-701-1033,
#' 01-701-1034, 01-701-1047, 01-701-1097, 01-705-1186,
#' 01-705-1292, 01-705-1310, 01-708-1286

USUBJID =
c("01-701-1015",
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
user = tibble(
  USUBJID = USUBJID)

result = inner_join(adlb,user)
admiral_adlb = result

# check
#USUBJID %in% result$USUBJID

# Finally, saved reduced ds
use_data(admiral_adlb, overwrite = TRUE)
