# R/admiral_verify_templates.R

# Create adam by sourcing template in admiral
source("inst/templates/ad_adae.R")

# Store adam data in cache.
theDir = tools::R_user_dir("admiral_templates_data", which = "cache")
load(paste0(theDir, "/adae.rda"))

# compare this adam to
library(pharmaverseadam )
adae2 = pharmaverseadam::adae
adae2

library(diffdf)

diffdf(adae, adae2 )
