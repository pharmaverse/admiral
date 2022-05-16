test_that("all templates are listed", {
  expect_identical(
    unclass(list_all_templates()),
    c("ADAE", "ADCM", "ADEG", "ADEX", "ADLB", "ADPP", "ADSL", "ADVS")
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

test_that("Error Message is returned if ADaM template file already exists", {
  dir <- tempdir()
  file <- file.path(dir, "adsl.R")
  use_ad_template("adsl", save_path = file, open = FALSE)

  expect_error(
    use_ad_template("adsl", save_path = file, open = FALSE)
  )
})

test_that("`convert_blanks_to_na.list` produces a lists", {
  x <- c("", "", "")
  expected_output <- lapply(x, convert_blanks_to_na)
  actual_output <- convert_blanks_to_na.list(x)

  expect_equal(expected_output, actual_output)
})
