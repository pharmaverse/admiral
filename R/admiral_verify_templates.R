# R/admiral_verify_templates.R

# TEMPLATES:   admiral/inst/templates/ad_{name}.R
# CACHE:     "/home/jim/.config/cache/R/admiral_templates_data

# Create adam by sourcing template in admiral
source("inst/templates/ad_adae.R")

# Store adam data in cache.
theDir = tools::R_user_dir("admiral_templates_data", which = "cache")

#  This makes all ds in pkg available
library(pharmaverseadam )
#
# need ?
# ~23 datasets in pharmaverseadam
original_datasets = data(package = "pharmaverseadam")
original_datasets = original_datasets$results
original_datasets = original_datasets[, "Title"]
original_datasets

adae2 = pharmaverseadam::adae
adcm2 = pharmaverseadam::adcm
`::`(pharmaverseadam, "adae")
`::`(pharmaverseadam, sym(original_datasets[[1]]))


do.call(`::`, list("pharmaverseadam", "adae"))
x="adae"
do.call(`::`, list("pharmaverseadam", x))


z=lapply(original_datasets, function(x)  do.call(`::`, list("pharmaverseadam", x)))
z

# find all templsates in `admiral` (when did admiral get loaded?)
templates = list.files("inst/templates", full.names=TRUE)
templates

# run the templates, store results tables in .cache (~23 tables)
lapply(templates, source)

# retrieve from cache
f = function(x = theDir){
unlist(lapply(list.files(theDir, full.names=TRUE), load))
}
theDir

ls()


original_datasets
templates


# compare ----
library(diffdf)


diffdf(adae, adae2 )
diffdf(adcm, adcm2 )
