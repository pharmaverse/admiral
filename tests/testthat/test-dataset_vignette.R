# dataset_vignette ----
## Test 1: A 'knitr_kable' object is outputted when run outside pkgdown ----
test_that("dataset_vignette Test 1: A 'knitr_kable' object is outputted when run outside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "false")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))

  dm <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~COUNTRY,
    "STUDY1", "1",      "USA",
    "STUDY1", "2",      "USA",
    "STUDY1", "3",      "USA",
    "STUDY1", "4",      "USA"
  )

  expect_s3_class(dataset_vignette(dm), "knitr_kable")
  expect_s3_class(dataset_vignette(dm, display_vars = exprs(STUDYID, USUBJID)), "knitr_kable")
})

## Test 2: A 'datatables' object is outputted when run inside pkgdown ----
test_that("dataset_vignette Test 2: A 'shiny.tag.list' is outputted when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))

  dm <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~COUNTRY,
    "STUDY1", "1",      "USA",
    "STUDY1", "2",      "USA",
    "STUDY1", "3",      "USA",
    "STUDY1", "4",      "USA"
  )


  expect_s3_class(dataset_vignette(dm), "shiny.tag.list")
  expect_s3_class(dataset_vignette(dm, display_vars = exprs(STUDYID, USUBJID)), "shiny.tag.list")
})

## Test 3: An error is outputted when calling variable not in dataset ----
test_that("dataset_vignette Test 3: An error is outputted when calling variable not in dataset", {
  dm <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~COUNTRY,
    "STUDY1", "1",      "USA",
    "STUDY1", "2",      "USA",
    "STUDY1", "3",      "USA",
    "STUDY1", "4",      "USA"
  )

  expect_error(dataset_vignette(dm, display_vars = exprs(AGE)))
})
