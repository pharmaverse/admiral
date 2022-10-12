library(admiral.test)

# dataset_vignette
## Test 1: A 'knitr_kable' object is output when run outside pkgdown ----
test_that("Test 1: A 'knitr_kable' object is output when run outside pkgdown", {
  expect_s3_class(dataset_vignette(head(admiral_dm)), "knitr_kable")
})

## Test 2: A 'datatables' object is output when run inside pkgdown ----
test_that("Test 2: A 'datatables' object is output when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(dataset_vignette(head(admiral_dm)), "datatables")
})

## Test 3: A 'knitr_kable' object is output when run outside pkgdown with display_vars ----
test_that("Test 3: A 'knitr_kable' object is output when run outside pkgdown with display_vars", {
  Sys.setenv(IN_PKGDOWN = "false")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(
    dataset_vignette(head(admiral_dm), display_vars = vars(STUDYID, USUBJID)),
    "knitr_kable"
  )
})

## Test 4: A 'datatables' object is output when run inside pkgdown with display_vars ----
test_that("Test 4: A 'datatables' object is output when run inside pkgdown with display_vars", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(
    dataset_vignette(head(admiral_dm), display_vars = vars(STUDYID, USUBJID)),
    "datatables"
  )
})
