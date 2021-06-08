context("test-convert_dtc_to_dtm")

test_that("Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_dtc_to_dtm("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:52")
  )
})
