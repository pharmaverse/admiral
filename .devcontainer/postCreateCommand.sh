#!/bin/bash

R -q -e 'renv::restore(); staged.dependencies::install_deps(staged.dependencies::dependency_table(project = ".", verbose = 1), verbose = 1);' >> post_create.log

jq --arg folder "$(pwd)/" '. + { "initial_working_directory": $folder }' .devcontainer/rstudio-prefs.json > ~/.config/rstudio/rstudio-prefs.json
