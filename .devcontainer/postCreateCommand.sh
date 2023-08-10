#!/bin/bash

R -q -e 'renv::restore(lockfile = file.path("renv", "profiles", paste(R.version$major, substr(R.version$minor, 1, 1), sep = "."), "renv.lock")); staged.dependencies::install_deps(staged.dependencies::dependency_table(project = ".", verbose = 1), verbose = 1);'

jq --arg folder "$(pwd)/" '. + { "initial_working_directory": $folder }' .devcontainer/rstudio-prefs.json > ~/.config/rstudio/rstudio-prefs.json
