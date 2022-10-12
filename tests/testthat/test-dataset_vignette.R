library(admiral.test)

# dataset_vignette
# ---- dataset_vignette, test 1: A 'knitr_kable' object is output when run outside pkgdown ----
test_that("dataset_vignette, test 1: A 'knitr_kable' object is output when run outside pkgdown", {
  expect_s3_class(dataset_vignette(head(admiral_dm)), "knitr_kable")
})

# ---- dataset_vignette, test 2: A 'datatables' object is output when run inside pkgdown ----
test_that("dataset_vignette, test 2: A 'datatables' object is output when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(dataset_vignette(head(admiral_dm)), "datatables")
})

# ---- dataset_vignette, test 3: A 'knitr_kable' object is output when run outside pkgdown using display_vars input ----
test_that("dataset_vignette, test 3: A 'knitr_kable' object is output when run outside pkgdown using display_vars input", {
  Sys.setenv(IN_PKGDOWN = "false")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(
    dataset_vignette(head(admiral_dm), display_vars = vars(STUDYID, USUBJID)),
    "knitr_kable"
  )
})

# ---- dataset_vignette, test 4: A 'datatables' object is output when run inside pkgdown using display_vars input ----
test_that("dataset_vignette, test 4: A 'datatables' object is output when run inside pkgdown using display_vars input", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(
    dataset_vignette(head(admiral_dm), display_vars = vars(STUDYID, USUBJID)),
    "datatables"
  )
})
