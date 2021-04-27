context("test-impute_dtc")

input <- c(
  "2019-07-18T15:25:40.243",
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

test_that("default: no date imputation, time part set o 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    as.character(NA),
    as.character(NA),
    as.character(NA)
  )
  expect_equal(impute_dtc(dtc = input), expected_output)
})

test_that("default: no date imputation,Missing time part imputed with 23:59:59 portion", {
  expected_output <- c(
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    as.character(NA),
    as.character(NA),
    as.character(NA)
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      time_imputation = "23:59:59"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc(
      dtc = input,
      time_imputation = "LAST"
    ),
    expected_output
  )
})

test_that("impute to first day/month if date is partial,Missing time part imputed with  00:00:00 portion", { # nolint
  expected_output <- c(
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-01T00:00:00",
    "2019-01-01T00:00:00",
    "2019-01-01T00:00:00"
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      date_imputation = "FIRST"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc(
      dtc = input,
      date_imputation = "01-01"
    ),
    expected_output
  )
})

test_that("impute to last day/month if date is partial,Missing time part imputed with 23:59:59 portion", { # nolint
  expected_output <- c(
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    "2019-02-28T23:59:59",
    "2019-12-31T23:59:59",
    "2019-12-31T23:59:59"
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      date_imputation = "LAST",
      time_imputation = "LAST"
    ),
    expected_output
  )
})

test_that("impute to MID day/month if date is partial,Missing time part imputed with 00:00:00 portion", { # nolint
  expected_output <- c(
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-15T00:00:00",
    "2019-06-15T00:00:00",
    "2019-06-15T00:00:00"
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      date_imputation = "MID"
    ),
    expected_output
  )
})
