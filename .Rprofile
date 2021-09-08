if (Sys.getenv("GITHUB_ACTIONS") == "" || (Sys.getenv("GITHUB_ACTIONS") == "true" && getRversion()$major == 3 && getRversion()$minor == 6)) {
  source("renv/activate.R")
} else {
  options(repos = c(CRAN = "https://cran.rstudio.com"))
}
