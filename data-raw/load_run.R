# PURPOSE:  gh actions uses this R code to setup/run script to run templates

# nolint start
# INSTALL all needed  packages
# load_all()
# SOURCE The script("admiral_verify_templates.R")
# RUN The function:  verify_templates()
# nolint end

install.packages("pak")
pak::pak("devtools")
pak::pak("diffdf")
# pak::pak("teal.data")


devtools::load_all()

source("data-raw/admiral_verify_templates.R") # nolint

# select templates to run
templates <- list_all_templates() |> paste()
templates <- templates[templates!="ADLBHY"]

# run only 2 templates
# verify_templates(ds = templates |> tolower() |> head(2))

# run all templates, default is all
verify_templates()
