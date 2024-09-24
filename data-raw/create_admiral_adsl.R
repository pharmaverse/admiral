#  Create dataset:   data/admiral_adsl.rda
#  This script:  create_admiral_adsl.R creates dataset data/admiral_adsl.rda.
#

# styler:style_pkg()
# devtools::lint()
# devtools::spellcheck()

# Preliminary
library(diffdf)


# To clarify directories (can be removed)
CACHE_DIR <- "~/.cache/R/admiral_templates_data/"
DATA_DIR <- "data-dir"
DATA_RAW <- "data-raw"
TEMPLATE_DIR <- "inst/templates/"

# clean CACHE_DIR
THE_FILE <- paste0(CACHE_DIR, "/adlb.rda")
THE_FILE <- paste0(CACHE_DIR, "/adsl.rda")
if (file.exists(THE_FILE)) file.remove(THE_FILE)

#
# STEPS
#
# First, use template to create the R script (in data-raw/admiral_adlb.R).
# Next, source this script and create the data (~/.cache/R/admiral_template_data/admiral_adlb.rda)
# Finally, shorten this data (now ~ 1.2 MB) by selecting only certain USUBJID

# OMIT -- orignal method  - OMIT
if (FALSE) {
  # ### original mehtod (method 1)
  # # First,  create the R script (from a template)
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
  load("~/.cache/R/admiral_templates_data/adlb.rda")
}

#
#
# Instead, USE template, as recommened by Buzz
#
source(paste0(TEMPLATE_DIR, "/ad_adsl.R")) # nolint
load(paste0(CACHE_DIR, "adsl.rda"))
admiral_adsl <- adsl

#
# Finally, save reduced ds
#
use_data(admiral_adsl, overwrite = TRUE)

#
#  TEST - is dataset identical to .... backup of unaltered dataset
#
e1 <- new.env()
e2 <- new.env()
load("data/admiral_adsl.rda", e1)

# CHANGE to YOUR location of original dataset
load("data-backup/admiral_adsl.rda", e2)


identical(e1$admiral_adsl, e2$admiral_adsl)
diffdf(compare=e1$admiral_adsl, base=e2$admiral_adsl, keys  = c("STUDYID", "USUBJID"))
capture.output(diffdf(compare=e1$admiral_adsl, base=e2$admiral_adsl, keys  = c("STUDYID", "USUBJID")),
               file="data-raw/diffdf_23SEPT")

#
# cleanup
#
#rm(e1)
#rm(e2)
