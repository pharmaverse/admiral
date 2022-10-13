# dataset_vignette ----
## Test 1: A 'knitr_kable' object is output when run outside pkgdown ----
test_that("dataset_vignette Test 1: A 'knitr_kable' object is output when run outside pkgdown", {
  dm <- tibble(
    STUDYID = rep("1001", 100),
    USUBJID = as.character(1:100),
    COUNTRY = rep("USA", 100)
  )

  expect_s3_class(dataset_vignette(dm), "knitr_kable")
})

## Test 2: A 'datatables' object is output when run inside pkgdown ----
test_that("dataset_vignette Test 2: A 'datatables' object is output when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))

  dm <- tibble(
    STUDYID = rep("1001", 100),
    USUBJID = as.character(1:100),
    COUNTRY = rep("USA", 100)
  )

  expect_s3_class(dataset_vignette(dm), "datatables")
})

## Test 3: A 'knitr_kable' object is output when run outside pkgdown with display_vars ----
test_that("dataset_vignette Test 3: A 'knitr_kable' object is output when run outside pkgdown with display_vars", {
  Sys.setenv(IN_PKGDOWN = "false")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))

  dm <- tibble(
    STUDYID = rep("1001", 100),
    USUBJID = as.character(1:100),
    COUNTRY = rep("USA", 100)
  )

  expect_s3_class(dataset_vignette(dm, display_vars = vars(STUDYID, USUBJID)), "knitr_kable")
})

## Test 4: A 'datatables' object is output when run inside pkgdown with display_vars ----
test_that("dataset_vignette Test 4: A 'datatables' object is output when run inside pkgdown with display_vars", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))

  dm <- tibble(
    STUDYID = rep("1001", 100),
    USUBJID = as.character(1:100),
    COUNTRY = rep("USA", 100)
  )

  expect_s3_class(dataset_vignette(dm, display_vars = vars(STUDYID, USUBJID)), "datatables")
})
