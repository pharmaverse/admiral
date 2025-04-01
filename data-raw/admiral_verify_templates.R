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
    get(ls()[ls() != "fileName"]) # get the new adlb dataset
}
save_rda <- function(data, file_path, new_name) {
  # new_name must include  .rda
  if (missing(new_name)) {
    save(data, file = file_path, compress="bzip2")
  }
}

compare = function(base, compare, keys, file){
   ##:ess-bp-start::conditional@		:##
browser(expr={		})##:ess-bp-end:##
  tryCatch({
      diffdf::diffdf(
        base = base,
        compare = compare,
        keys = keys,
        file = paste0(dataset_dir,"/diff_", file, ".txt"))
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
rlang::global_entrace()
clean_cache()

library(pharmaverseadam)

# SOURCE:  Use pharmaverseadam as source for ADaM *.rda files (22 ADaM files) 
# ignore adlbhy.rda (per Ben) 


# REMOVE loop:   per Ben, only do "admiral"
packages_list = c("admiral")
for (pkg in packages_list){                  # ---- pkg ----
  sprintf("generating ADaMs for  %s package", pkg)
  library(pkg, character.only = TRUE)

  # DISCUSS: Ask user to update <pkg>?

  # per Ben, ignore *.rda files in admiral/data
  source_adams=data(package= "pharmaverseadam")

  # 23 ADaMs found
  source_adams = source_adams$results[,"Item"]

  # per Ben, ignore "adlbhy"
  #source_adams = source_adams[source_adams != "adlbhy"]
  source_adams

  # gather al templates for this pkg (12 found)
  template_path = file.path(system.file(package = pkg), "templates")
  templates = list.files(template_path, pattern = "ad_")
  templates

  ## BUT, Fewer templates than  ADaMs in pharmaverseadam!
  ## Reduce number of source_adams ! 

  
  ## templates can generate these ADaMs, 12
  names = sapply(templates, function(x) gsub("ad_|\\.R","",x ), USE.NAMES = FALSE)
  ## name our templates



  # KEEP only elments in both (now 12)
source_adams =  source_adams[source_adams %in% names]
  # finally name the templates vector
  names(templates) = names 
  templates


  # new, generated ADaM datasets will be put in cache dir
#  dataset_dir = tools::R_user_dir(sprintf("%s_templates_data", pkg), which="cache")
#  dataset_dir = tools::R_user_dir(sprintf("%s_templates_data", pkg), which="cache")

  # collect in list of important paths
  path = list("templates_path" = file.path(system.file(package = pkg), "templates"),
              dataset_dir = tools::R_user_dir("admiral_templates_data", which="cache") 
  )

  ## TODO:  TESTING:   ignore all templates but these two 
  ignore_templates_pkg = templates[!(templates %in% c("ad_adlb.R", "ad_adsl.R"))]

  ## AFTER Testing, this is only template to be ignored
  #ignore_templates_pkg = templates[!(templates %in% c("ad_adlbhy.R"))]

    # begin tp  ----
    for (tp in templates) {
      if (tp == "ad_adlb.R") next  ## adlb is big, skip while debugging
      # Each tp creates single ADaM package

      # 
        name = gsub("ad_|\\.R", "", tp)
        adam_old = do.call(`::`, args=list("pharmaverseadam", name))
      # run template, which caches new adam_new
      if(tp %in% ignore_templates_pkg) {
        cat("Ignoring template: ", tp, "\n")
        next
      } else {
        cat("Running template: ", tp, "\n")
        source(file.path(template_path, tp), echo = TRUE) # nolint
      }
      # retrieve new adam from cache dir (puts into globalenv())
      load_rda(paste0(dataset_dir,"/", name, ".rda"))
      browser()
      adam_new = get(name)

      # DISCUSS
      all_results <- c()

##   DISCUSS:   Ben:  use full dataset, no reduction

       ## DISCUSSION:   adlb.rda ADaM is large and must be reduced in size (for testing)
       ## What about other ADaMs ?  use vector
      if (FALSE) {
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

       }
          # Finally, save  dataset (reduced or not) in data/
          # TODO:  use vector
          if (tp == "ad_adlb.R") {
            save_rda(adam_new, file_path="data/admiral_adlb.rda")
          } else {
            #save_rda(adam_new, file_path="data/admiral_adsl.rda")
            save_rda(adam_new, file_path="~/admiral_adsl.rda")
          }
       } ## end if

    ## TODO function use vectors
  browser()
     sprintf("comparing ... diffdf")
    #sprintf("comparing ....%s", tp)
    cat("comparing ...", tp, "\n")
    compare(base=adam_old,
             compare=adam_new,
            keys= ifelse(tp == "adlb.R",
                         c("USUBJID", "PARAMCD", "AVISIT", "ADT"),
                         c("STUDYID","USUBJID"),
             file = paste0("data/diff_", name, ".txt")
                         )

} ## end tp
} ## end for packages_list

print("DONE")
