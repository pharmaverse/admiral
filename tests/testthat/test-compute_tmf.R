context("test-compute_tmf")

test_that("compute TMF", {
  input_dtc <- c(
    "2019-07-18T15:25:52",
    "2019-07-18T15:25",
    "2019-07-18T15",
    "2019-07-18",
    "2019-02",
    "2019",
    "2019---07"
  )
  input_dtm <- c(
    as.POSIXct("2019-07-18T15:25:52"),
    as.POSIXct("2019-07-18T15:25:00"),
    as.POSIXct("2019-07-18T15:00:00"),
    as.POSIXct("2019-07-18"),
    as.POSIXct("2019-02-01"),
    as.POSIXct("2019-01-01"),
    as.POSIXct("2019-01-01")
  )
  expected_output <- c(
    "",
    "S",
    "M",
    "H",
    "H",
    "H",
    "H"
  )

  expect_equal(
    compute_tmf(
      dtc = input_dtc,
      dtm = input_dtm
    ),
    expected_output
  )
})
