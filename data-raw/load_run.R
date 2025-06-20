# This script:  data-raw/load_all.R
# PURPOSE:  gh actions uses this R code.

# nolint start
# install all needed  packages
# load_all()
# source("admiral_verify_templates.R")
# run the function:  verify_templates()
# nolint end

install.packages("pak")
pak::pak("devtools")
pak::pak("diffdf")
pak::pak("teal.data")


devtools::load_all()

source("data-raw/admiral_verify_templates.R") # nolint

templates <- list_all_templates() |> paste()
templates <- templates[templates!="ADLBHY"]

verify_templates(ds = templates |> tolower() |> head(2))
