context("test-convert_date_to_dtm")

test_that("Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18T15:25:52"),
    as_iso_dtm("2019-07-18T15:25:52")
  )
})

test_that("Impute incomplete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18", time_imputation = "23:59:59"),
    as_iso_dtm("2019-07-18T23:59:59")
  )
})

test_that("Convert -- DT into a date time object", {
  expect_equal(
    convert_date_to_dtm(as.Date("2019-07-18"), time_imputation = "23:59:59"),
    as_iso_dtm("2019-07-18T23:59:59")
  )
})

test_that("Keep -- DTM as the original date time object", {
  expect_equal(
    convert_date_to_dtm(as_iso_dtm("2019-07-18T15:25:52"), time_imputation = "23:59:59"),
    as_iso_dtm("2019-07-18T15:25:52")
  )
})
