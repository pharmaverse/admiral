# This script:  data-raw/load_all.R

#
# install all needed  packages
# load_all()
# source("admiral_verify_templates.R")
# run the function:  verify_templates():w

install.packages("pak")
pak::pak("devtools")
pak::pak("diffdf")
pak::pak("teal.data")


devtools::load_all()

source("data-raw/admiral_verify_templates.R")

verify_templates()










