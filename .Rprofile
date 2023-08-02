# Set renv profile base on R version.
renv_profile <- paste(R.version$major, substr(R.version$minor, 1, 1), sep = ".")
if (file.exists("./renv/profiles")) {
  message("Using renv profile from `renv/profile` file.")
  if (renv_profile %in% c("4.1", "4.2", "4.3")) {
    message("Set renv profile to `", renv_profile, "`")
    Sys.setenv("RENV_PROFILE" = renv_profile)
  } else {
  message("This repository do not contains the renv profile for your R version.")
  }
}

if ((Sys.getenv("GITHUB_ACTIONS") != "") || (Sys.getenv("DOCKER_CONTAINER_CONTEXT") != "")) {
  options(repos = c(CRAN = "https://cran.rstudio.com"))
  Sys.setenv(RENV_AUTOLOADER_ENABLED=FALSE)
}
source("renv/activate.R")