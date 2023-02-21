file <- "/cloud/project/inst/templates/ad_admh.R"
#
templates <- list.files(
  "/cloud/project/inst/templates/",
  pattern = "\\.R$",
  ignore.case = TRUE,
  full.names = TRUE
)

str = gsub("([\\])","", str)

cmd <- sprintf("Rscript --vanilla %s", file)
exe <- "Rscript --vanilla -e"
wchl <- "'withCallingHandlers(source('"
wchr <- "\\\\'), lifecycle_warning_deprecated = function(e){stop(e)})'"
cmd <- stringr::str_c(exe, " ", wchl, file, wchr)

system(cmd)

s1 <- "///\\\"

Rscript --vanilla -e "withCallingHandlers(source(\"/cloud/project/inst/templates/ad_admh.R\"), lifecycle_warning_deprecated = function(e){stop(e)})"
#Rscript --vanilla -e withCallingHandlers(source("/cloud/project/inst/templates/ad_admh.R\"), lifecycle_warning_deprecated = function(e){stop(e)})
