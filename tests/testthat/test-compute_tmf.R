context("test-compute_tmf")

inputdtc <- c(
  "2019-07-18T15:25:52",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)
inputdtm <- c(
  is.POSIXct("2019-07-18T15:25:52"),
  is.POSIXct("2019-07-18T15:25:00"),
  is.POSIXct("2019-07-18T15:00:00"),
  is.POSIXct("2019-07-18"),
  is.POSIXct("2019-02-01"),
  is.POSIXct("2019-01-01"),
  is.POSIXct("2019-01-01")
)
test_that("compute TMF", {
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
      dtc = inputdtc,
      dtm = inputdtm
    ),
    expected_output
  )
})
