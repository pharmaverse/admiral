
# -------   N ----
# Step 1, gather 23 datasets from pharamaversadam ---

#  all ds available ---
library(pharmaverseadam )
# ...

# 23- get datasets from pharmaverseadam package ("original")
old = data(package = "pharmaverseadam")
old = old$results
old = old[, "Title"]
old

# but only use 6 ?
old = list(adae, adcm, adeg, adex, adlb, adsl)
old[[1]]


# step 2 - generate new datasets with admiral templates  ----

# 12-find all templates in `admiral` (when did admiral get loaded?)
templates = list.files("inst/templates", full.names=TRUE)
templates

# But, instead of 12, only want 6

files = c("adae.R", "adcm.R", "adeg.R", "adex.R", "adlb.R", "adsl.R")
templates = paste0("inst/templates/ad_", files)
templates

# 6, slow ...
lapply(templates, source)

# check .rda from cache, only 6
theDir = tools::R_user_dir("admiral_templates_data", which = "cache")
theDir
list.files(theDir)

# What is in the cache?
new = paste0(theDir,"/", dir(theDir))
new

# load new files from  cache, but in separate environment
newE = new.env()
lapply(new, load, newE)
ls(newE)
newE$adlb

# compare ----
library(diffdf)
diffdf(
  base = old,
  compare = adae, # generated
  keys = c("")
)

old[[1]]
new = new.env()
list.files(theDir)[[1]]
load(paste0(theDir,"/", "adae.rda"))
ls(new)
new$adae

diffdf(
  base = old[[1]],   # adae
  compare =newE$adae
)

diffdf(
  base = old[[2]],
  compare = newE$adcm
)



