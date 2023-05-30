# Set renv profile base on R version.
renv_profile <- paste(R.version$major, substr(R.version$minor, 1,1), sep = ".")
if (file.exists("./renv/profile")){
  message("Using renv profile from `renv/profile` file.")
} else if (renv_profile %in% c("4.1", "4.2", "4.3")) {
  message("Set renv profile to", renv_profile)
  Sys.setenv("RENV_PROFILE" = renv_profile)
} else {
  message("Using default renv profile")
}

if (Sys.getenv("GITHUB_ACTIONS") == "" || (Sys.getenv("GITHUB_ACTIONS") == "true" && getRversion()$major == 3 && getRversion()$minor == 6)) {
  source("renv/activate.R")
} else {
  options(repos = c(CRAN = "https://cran.rstudio.com"))
}


