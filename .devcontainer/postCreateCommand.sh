#!/bin/bash

R -q -e 'renv::activate(profile = paste(R.version$major, substr(R.version$minor, 1, 1), sep = ".")); renv::restore(); staged.dependencies::install_deps(staged.dependencies::dependency_table(project = ".", verbose = 1), verbose = 1);'

jq --arg folder "$(pwd)/" '. + { "initial_working_directory": $folder }' .devcontainer/rstudio-prefs.json > ~/.config/rstudio/rstudio-prefs.json
