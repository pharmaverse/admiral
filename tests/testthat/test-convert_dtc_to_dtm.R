context("test-convert_dtc_to_dtm")

inputdtc <- c(
  "2019-07-18T15:25:52"
)


test_that("Convert a complete -- DTC into a date time object", {
  expected_output <- c(
    ymd_hms("2019-07-18T15:25:52")
  )
  expect_equal(
    convert_dtc_to_dtm(dtc = inputdtc),
    expected_output
  )
})
