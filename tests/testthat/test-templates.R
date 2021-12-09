test_that("template scripts execute without errors", {
  old_wd <- getwd()
  setwd("../..")
  on.exit(setwd(old_wd))

  templates <- list.files(
    system.file("templates", package = "admiral"),
    full.names = TRUE
  )

  for (file in templates) {
    expect_error(
      suppressWarnings(source(file, local = new.env())),
      NA,
      info = sprintf("Template script %s", basename(file))
    )
  }
})
