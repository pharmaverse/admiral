#
# Test:
#    In folder data/   files should contain identical R objects:

#    example_qs.rda was created in via script inst/example_scripts/example_qs.R
#    example_qs_new.rda was created by NEWER method, script data-raw/example_qs_new.R

# Run both methods, create the 2 *.rda files:
source("data-raw/example_qs.r")
source("inst/example_scripts/example_qs.R")


## Please restart R or remove objects in environment

# load stored objects into environment
# original method
load("data/example_qs.rda")

# newer method
load("data/example_qs_new.rda")


identical(example_qs_new, example_qs)
