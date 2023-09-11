# Set renv profile base on R version.
.get_dependencies <- function(project_dir) {

  admdev_loc <- find.package("admiraldev", lib.loc = .libPaths(), quiet = TRUE)
  adm_dev_suggests <- if(length(admdev_loc) != 0) {
    renv:::renv_dependencies_discover_description(admdev_loc, fields = c("Depends", "Imports", "LinkingTo", "Suggests"))
  } else {
    data.frame(Packages = character(0))
  }
  suggests_packages <- renv:::renv_dependencies_discover_description(project_dir, fields = "Suggests")

  packages <- names(
    renv:::renv_package_dependencies(
      unique(c(
        project_dir,
        adm_dev_suggests[["Package"]],
        suggests_packages[["Package"]],
        c("staged.dependencies", "renv", "styler")
      ))
    )
  )
  packages[!(packages %in% c("admiral", "admiraldev", "admiralci", "admiral.test", "pharmaversesdtm", getwd()))]
}

options(renv.snapshot.filter = .get_dependencies)

.renv_profile <- paste(R.version$major, substr(R.version$minor, 1, 1), sep = ".")
if (!file.exists("./renv/profile")) {
  if (.renv_profile %in% c("4.1", "4.2", "4.3")) {
    message("Set renv profile to `", .renv_profile, "`")
    Sys.setenv("RENV_PROFILE" = .renv_profile)
  } else {
    message("This repository do not contains the renv profile for your R version.")
  }
} else {
  message(
    "Using renv profile from `renv/profile` file.\n",
    "The `", readLines("./renv/profile"), "` profile will be used."
  )
}

if (Sys.getenv("GITHUB_ACTIONS") != "") {
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
  Sys.setenv("RENV_AUTOLOADER_ENABLED" = FALSE)
}
Sys.setenv("RENV_CONFIG_SANDBOX_ENABLED" = FALSE)
Sys.setenv("RENV_CONFIG_AUTO_SNAPSHOT" = FALSE)
source("renv/activate.R")
