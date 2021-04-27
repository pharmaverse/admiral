context("test-compute_dtf")

inputdtc <- c(
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)
inputdt <- c(
  as.Date("2019-07-18"),
  as.Date("2019-02-01"),
  as.Date("2019-01-01"),
  as.Date("2019-01-01")
)
test_that("compute DTF", {
  expected_output <- c(
    NA,
    "D",
    "M",
    "M"
  )
  expect_equal(
    compute_dtf(
      dtc = inputdtc,
      dt = inputdt
    ),
    expected_output
  )
})
