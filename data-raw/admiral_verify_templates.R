#  This script:  data-raw/admiral_verify_templates.R
#  Create ADaM datasets:   data/<package>_<name>.rda (ex: data/admiral_adlb.rda)

# <package> = name of package (ex:  "admiral )
# <name> = ADaM name  (ex:  adlb)
# <name>_old =  variable for original ADaM dataset before we anything (ex: adlb_old)
# ad_<name>.R  template names, where name is 4 more characters (ex: ad_adlb.R)

# data/<package>_<name>.rda as it in admiral/data before do anything (ex: data/admiral_adlb.rda)
# adlb      after template is run, before reduced
# admiral_adlb after template is run, reduced and saved


# LOCATIONS
# dataset_dir (cache)     ~/.config/R/admiral_templates_data   (varies by OS/configuration)
# admiral/data     dir where admiral_adlb.rda will be saved
# admiral/inst/templates/ template to create ADaM datasets
#

#

## packages_list
packages_list = c("admiral")

for (pkg in packages_list){

  # 1st, make copy of all ADaMs in pkg (TODO:  all ADaMs)
  adlb_old <- admiral::admiral_adlb

  # setup
  sprintf("generating ADaMs for  %s package", pkg)
  assign(x=pkg, value = pkg)
  library(pkg, character.only = TRUE)


  # update <pkg>?    ASSUME:  already done

  # gather all package templates
  template_path = file.path(system.file(package = pkg), "templates")
  templates = list.files(template_path, pattern = "ad_")

  # TODO:  Construct named character vector, one name for each pkg,  each referring to have 0,1  or more templates to ignore
  # c("admiral" = c("ad_adlb.R", "ad_adlb2.R", ...), "package2" = c("ad_adlb.R", "ad_adlb2.R", ...)))

  ## TODO:   templates to ignore, fill in MANUALLY
  ignore_templates= c("admiral" = c())

  ## TODO:  for this pkg, which templates to ignore
  #ignore_templates_pkg = ignore_templates[ignore_templates == pkg)

  ## TODO:   for now, use this, just one template, ignore the others
  ignore_templates_pkg = templates[templates != "ad_adlb.R"]

    for (tp in templates) {
      if(tp %in% ignore_templates_pkg) {
        cat("Ignoring template: ", tp, "\n")
      } else {
        cat("Running template: ", tp, "\n")
        #source(file.path(template_path, tp), echo = TRUE) # nolint
      }
    }  ## end for

   # new ADaM datasets are now in this cache dir
  (dataset_dir = tools::R_user_dir(sprintf("%s_templates_data", pkg), which="cache"))
  dir(dataset_dir)

  ## TODO: load all the new ADaM datasets
  # for now, put into globalenv() as admiral_adlb
  load(paste0(dataset_dir, "/adlb.rda"))

}  ## end for all packages


all_results <- c()

## DISCUSSION:   adlb.rda ADaM is large and must be reduced in size (for testing)
## What about other ADaMs ?

# Limit rows by selecting only these USUBJIDs
usubjids <-
  c(
    "01-701-1015",
    "01-701-1023",
    "01-701-1028",
    "01-701-1033",
    "01-701-1034",
    "01-701-1047",
    "01-701-1097",
    "01-705-1186",
    "01-705-1292",
    "01-705-1310",
    "01-708-1286"
  )

# small adlb
admiral_adlb <- dplyr::filter(adlb, USUBJID %in% usubjids)



# Finally, save reduced dataset in data/
usethis::use_data(admiral_adlb, overwrite = TRUE)

## TODO function to go compare all the new ADaM against the old one.
## TODO:  use safely to catch errors ?    Or tryCatch?
# Compare with previous version

diffdf::diffdf(
  base = adlb_old,
  compare = admiral_adlb,
  keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT")
)

