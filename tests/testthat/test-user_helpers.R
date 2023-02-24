# list_all_templates ----
test_that("all templates are listed", {
  expect_equal(
    unclass(list_all_templates()),
    c("ADAE", "ADCM", "ADEG", "ADEX", "ADLB", "ADLBHY", "ADMH", "ADPC", "ADPP", "ADSL", "ADVS"),
    ignore_attr = TRUE
  )
})

test_that("Error Message is returned if package is not installed", {
  expect_error(
    list_all_templates(package = "non-existing-package"),
    regexp = paste0(
      "No package called 'non-existing-package' is installed ",
      "and hence no templates are available"
    )
  )
})

# use_ad_template ----
test_that("package templates can be used", {
  dir <- tempdir()
  suppressMessages(file <- file.path(dir, "advs.R"))
  suppressMessages(use_ad_template("advs", save_path = file, open = FALSE))

  expect_true(file.exists(file))
  expect_identical(
    readLines(system.file("templates/ad_advs.R", package = "admiral")),
    readLines(file)
  )
  file.remove(file)
})

test_that("Error Message is returned if no ADaM template is available", {
  dir <- tempdir()
  suppressMessages(file <- file.path(dir, "adxx.R"))
  expect_error(
    suppressMessages(use_ad_template("adxx", save_path = file, open = FALSE)),
    "No template for 'ADXX' available."
  )
})

test_that("Error Message is returned if ADaM template file already exists", {
  dir <- tempdir()
  suppressMessages(file <- file.path(dir, "adsl.R"))
  suppressMessages(use_ad_template("adsl", save_path = file, open = FALSE))

  expect_error(
    suppressMessages(use_ad_template("adsl", save_path = file, open = FALSE))
  )
  file.remove(file)
})

# print.adam_templates ----
test_that("`adam_templates` objects are printed as intended: no templates", {
  templates <- list_all_templates(package = "dplyr")
  expected_print_output <- c(
    "No ADaM templates available in package 'dplyr'"
  )
  expect_identical(capture.output(print(templates)), expected_print_output)
})

test_that("`adam_templates` objects are printed as intended: some templates", {
  templates <- c("ADAE", "ADSL") %>%
    structure(class = c("adam_templates", "character"), package = "admiral")
  expected_print_output <- c(
    "Existing ADaM templates in package 'admiral':",
    if (is.na(iconv("\U2022"))) "- ADAE" else "\U2022 ADAE",
    if (is.na(iconv("\U2022"))) "- ADSL" else "\U2022 ADSL"
  )
  expect_identical(capture.output(print(templates)), expected_print_output)
})
