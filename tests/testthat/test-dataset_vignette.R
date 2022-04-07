test_that("a 'knitr_kable' object is output when run outside pkgdown", {
  expect_s3_class(dataset_vignette(head(adsl)), "knitr_kable")
})

test_that("a 'datatables' object is output when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(dataset_vignette(head(adsl)), "datatables")
})

test_that("a 'knitr_kable' object is output when run outside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "false")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(
    dataset_vignette(head(adsl), display_vars = vars(STUDYID, USUBJID)),
    "knitr_kable")
})
