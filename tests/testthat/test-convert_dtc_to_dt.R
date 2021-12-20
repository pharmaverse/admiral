inputdtc <- c(
  "2019-07-18T15:25:52",
  "2019-07-18"
)

test_that("Convert a complete -- DTC into a date object", {
  expected_output <- c(
    as.Date("2019-07-18"),
    as.Date("2019-07-18")
  )
  expect_equal(
    convert_dtc_to_dt(dtc = inputdtc),
    expected_output
  )
})
