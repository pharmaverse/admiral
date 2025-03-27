#  This script:  data-raw/admiral_verify_templates.R
#  Create ADaM datasets:   data/<package>_<name>.rda (ex: data/admiral_adlb.rda)

#' @parma pkg  package (ex:  "admiral )
#' @param name ADaM name, without prefix or suffix  (ex:  adlb)
#' @param adam_new  ADaM after template is run, reduced and saved
#' @param adam_old =  original ADaM dataset before we anything

#' @param tp   template R file  (ex: ad_adlb.R)
#' @templates_path path to templates where package is installed (not sourced)
#' admiral/inst/templates/ template to create ADaM datasets
#' @param dataset_dir (cache)     ~/.config/R/admiral_templates_data   (varies by OS/configuration)
# admiral/data     dir where admiral_adlb.rda will be saved

# cache (1st version)
load_rda = function(fileName) {
    load(fileName)
    get(ls()[ls() != "adlb"]) # get the new adlb dataset
}
# better version
load_rda2 <- function (path,  name){
 # in cache dir (adlb.rda)
  load(paste0(path,"/", name, ".rda"), envir=globalenv(), verbose=TRUE)
}


save_rda = function(data, file_path, new_name) {
  # new_name must include  .rda
  if (missing(new_name)) {
    save(data, file = file_path, compress="bzip2")
  }
}

compare = function(base, compare, keys, file){
    tryCatch({
      diffdf::diffdf(
        base = base,
        compare = compare,
        keys = keys,
        file = paste0(dataset_dir,"/diff_", name, ".txt")
     )
  },
      error = function(e) {
        message("Error in diffdf: ", e$message)
  }
)  ## end tryCatch
}  ## end compare

clean_cache = function() {
  dataset_dir = tools::R_user_dir("admiral_templates_data", which="cache")
  if (dir.exists(dataset_dir)) {
    unlink(dataset_dir, recursive = TRUE)
    message("Cache directory deleted: ", dataset_dir)
  } else {
    message("Cache directory does not exist: ", dataset_dir)
  }
}

main = function() {

clean_cache()
packages_list = c("admiral")

for (pkg in packages_list){                  # ----
  sprintf("generating ADaMs for  %s package", pkg)
  library(pkg, character.only = TRUE)

  # DISCUSS: ask user to update <pkg>?

  # gather all package templates
  template_path = file.path(system.file(package = pkg), "templates")
  templates = list.files(template_path, pattern = "ad_")

  # new ADaM datasets will be put in cache dir
   dataset_dir = tools::R_user_dir(sprintf("%s_templates_data", pkg), which="cache")

  ## DISCUSS:   templates to ignore, add  MANUALLY?
  ## TODO:   for now, ignore templates EXCEPT these 2
  ignore_templates_pkg = templates[!(templates %in% c("ad_adlb.R", "ad_adsl.R"))]

    # begin tp  ----
    for (tp in templates) {
      #if (tp == "ad_adlb.R") next  ## adlb is big, skip while debugging
      # Each tp creates single ADaM package

      # 1st, make copy of all ADaMs in pkg
      # DISCUSS:  put in separate environment?

      # (for now) use:   adam_old
      # TODO:  automate, use vector

      if (tp == "ad_adlb.R") {
        adam_old <- admiral::admiral_adlb
      } else {
        adam_old <- admiral::admiral_adsl
      }

      # extract base ADaM name (ex:  adsl)
      # TODO: use vector
      name = gsub("ad_", "", tp)
      name = gsub(".R", "", name)
      name

      # run template, which caches new adam_new
      if(tp %in% ignore_templates_pkg) {
        cat("Ignoring template: ", tp, "\n")
        next
      } else {
        cat("Running template: ", tp, "\n")
        source(file.path(template_path, tp), echo = TRUE) # nolint
      }
      # retrieve new adam from cache dir (puts into globalenv())
      load_rda2(dataset_dir, name)

      all_results <- c()

       ## DISCUSSION:   adlb.rda ADaM is large and must be reduced in size (for testing)
       ## What about other ADaMs ?  use vector
       if (tp == "ad_adlb.R") {
         # Reduce size of adlb dataset for testing
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
            adam_new = admiral_adlb
       } else {
             admiral_adsl <- adsl
             adam_new = admiral_adsl
       }  ## end if

          # Finally, save  dataset (reduced or not) in data/
          # TODO:  use vector
          if (tp == "ad_adlb.R") {
            save_rda(adam_new, file_path="data/admiral_adlb.rda")
          } else {
            save_rda(adam_new, file_path="data/admiral_adsl.rda")
          }

    ## TODO function use vectors
    sprintf("comparing ....%s", tp)
    cat("comparing ...", tp, "\n")
    compare(base=adam_old,
             compare=adam_new,
             keys= ifelse(tp == "adlb.R", c("USUBJID", "PARAMCD", "AVISIT", "ADT"), NULL),
             file = paste0("data/diff_", name, ".txt")
                         )

} ## end tp
} ## end for packages_list

print("DONE")
}
