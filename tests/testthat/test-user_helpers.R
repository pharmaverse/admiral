# list_all_templates ----
## Test 1: all templates are listed ----
test_that("list_all_templates Test 1: all templates are listed", {
  expect_equal(
    unclass(list_all_templates()),
    c(
      "ADAE", "ADCM", "ADEG", "ADEX", "ADLB", "ADLBHY", "ADMH", "ADPC", "ADPP",
      "ADPPK", "ADSL", "ADVS"
    ),
    ignore_attr = TRUE
  )
})

## Test 2: Error Message is returned if package is not installed ----
test_that("list_all_templates Test 2: Error Message is returned if package is not installed", {
  expect_snapshot(
    list_all_templates(package = "non-existing-package"),
    error = TRUE
  )
})

# use_ad_template ----
## Test 3: package templates can be used ----
test_that("use_ad_template Test 3: package templates can be used", {
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

## Test 4: Error Message is returned if no ADaM template is available ----
test_that("use_ad_template Test 4: Error Message is returned if no ADaM template is available", {
  dir <- tempdir()
  suppressMessages(file <- file.path(dir, "adxx.R"))
  expect_snapshot(
    suppressMessages(use_ad_template("adxx", save_path = file, open = FALSE)),
    error = TRUE
  )
})

## Test 5: Error Message is returned if ADaM template file already exists ----
test_that("use_ad_template Test 5: Error Message is returned if ADaM template file already exists", {
  dir <- tempdir()
  suppressMessages(file <- file.path(dir, "adsl.R"))
  suppressMessages(use_ad_template("adsl", save_path = file, open = FALSE))

  expect_snapshot(
    suppressMessages(use_ad_template("adsl", save_path = file, open = FALSE)),
    error = TRUE
  )
  file.remove(file)
})

# print.adam_templates ----
## Test 6: `adam_templates` objects are printed as intended: no templates ----
test_that("print.adam_templates Test 6: `adam_templates` objects are printed as intended: no templates", {
  templates <- list_all_templates(package = "dplyr")
  expected_print_output <- c(
    "No ADaM templates available in package 'dplyr'"
  )
  expect_identical(capture.output(print(templates)), expected_print_output)
})

## Test 7: `adam_templates` objects are printed as intended: some templates ----
test_that("print.adam_templates Test 7: `adam_templates` objects are printed as intended: some templates", {
  templates <- c("ADAE", "ADSL") %>%
    structure(class = c("adam_templates", "character"), package = "admiral") # nolint: undesirable_function_linter
  expected_print_output <- c(
    "Existing ADaM templates in package 'admiral':",
    if (is.na(iconv("\U2022"))) "- ADAE" else "\U2022 ADAE",
    if (is.na(iconv("\U2022"))) "- ADSL" else "\U2022 ADSL"
  )
  expect_identical(capture.output(print(templates)), expected_print_output)
})
