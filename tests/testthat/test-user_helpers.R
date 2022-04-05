test_that("all templates are listed", {
  expect_identical(
    unclass(list_all_templates()),
    c("ADAE", "ADCM", "ADEG", "ADEX", "ADSL", "ADVS")
  )
})

test_that("package templates can be used", {
  dir <- tempdir()
  file <- file.path(dir, "advs.R")
  use_ad_template("advs", save_path = file, open = FALSE)

  expect_true(file.exists(file))
  expect_identical(
    readLines(system.file("templates/ad_advs.R", package = "admiral")),
    readLines(file)
  )
})

test_that("Error Message is returned if no ADaM template is available", {
  dir <- tempdir()
  file <- file.path(dir, "adxx.R")
  expect_error(
    use_ad_template("adxx", save_path = file, open = FALSE),
    "No template for 'ADXX' available."
  )
})

