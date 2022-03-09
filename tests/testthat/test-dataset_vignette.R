test_that("a 'knitr_kable' object is output when run outside pkgdown", {
  expect_s3_class(dataset_vignette(head(adsl)), "knitr_kable")
})

test_that("a 'datatables' object is output when run inside pkgdown", {
  Sys.setenv(IN_PKGDOWN = "true")
  on.exit(Sys.setenv(IN_PKGDOWN = ""))
  expect_s3_class(dataset_vignette(head(adsl)), "datatables")
})
