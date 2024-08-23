
#
# Test:
#    Newly created *.rda files in folder data/   MUST BE IDENTICAL to old_data/*.rda

#    *.rda file in old_data/ is copy of file created in via script inst/example_scripts/example_qs.R

# Run both methods, create the 2 *.rda files:
source("data-raw/example_qs.R")  # new
source("inst/example_scripts/example_qs.R")  # old


## Please restart R or remove objects in environment

# Load stored objects into environments
# Because `load` places objects in current environment, put old_data in a new environment

# For created data, place object into current_env
load("data/example_qs.rda")

# But create new environment for old_data/
e <- new_environment()
load("old_data/example_qs.rda", envir = e)


# test
identical(example_qs, e$example_qs)

# cleanup: remove e
rm(e)
