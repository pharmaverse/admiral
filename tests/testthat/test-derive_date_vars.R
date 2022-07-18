library(lubridate)

input <- c(
  "2019-07-18T15:25:40.243",
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07",
  "2003-12-15T-:15:18",
  "2003-12-15T13:-:19",
  "2020-07--T00:00"
)

# impute_dtc ----
## Test 1: default: no date imputation, time part set to 00:00:00 ----
test_that("impute_dtc Test 1: default: no date imputation, time part set to 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    NA_character_,
    NA_character_,
    NA_character_,
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    NA_character_
  )
  expect_equal(impute_dtc(dtc = input), expected_output)
})

## Test 2: no date imputation, min and sec imputed with 59 ----
test_that("impute_dtc Test 2: no date imputation, min and sec imputed with 59", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "2003-12-15T13:59:59",
    NA_character_
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "23:59:59"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "LAST"
    ),
    expected_output
  )
})

## Test 3: impute date to first, time to 00:00:00 ----
test_that("impute_dtc Test 3: impute date to first, time to 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-01T00:00:00",
    "2019-01-01T00:00:00",
    "2019-01-01T00:00:00",
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    "2020-07-01T00:00:00"
  )

  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "01-01"
    ),
    expected_output
  )
})

## Test 4: impute day to last, time to 23:59:59 and preserve = FALSE ----
test_that("impute_dtc Test 4: impute day to last, time to 23:59:59 and preserve = FALSE", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    "2019-02-28T23:59:59",
    NA_character_,
    NA_character_,
    "2003-12-15T23:59:59",
    "2003-12-15T13:59:59",
    "2020-07-31T23:59:59"
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "D",
      date_imputation = "LAST",
      time_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 5: impute date to last, time to 23:59:59 and preserve = TRUE ----
test_that("impute_dtc Test 5: impute date to last, time to 23:59:59 and preserve = TRUE", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    "2019-02-28T23:59:59",
    "2019-12-31T23:59:59",
    "2019-12-07T23:59:59",
    "2003-12-15T23:15:18",
    "2003-12-15T13:59:19",
    "2020-07-31T00:00:59"
  )
  expect_equal(
    imputes <- impute_dtc(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "LAST",
      time_imputation = "LAST",
      preserve = TRUE
    ),
    expected_output
  )
})

## Test 6: impute date to last, time to 23:59:59 and preserve = FALSE ----
test_that("impute_dtc Test 6: no date imputation, impute second 59 and preserve = FALSE", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  )

  expect_equal(
    imputes <- impute_dtc(
      dtc = input,
      highest_imputation = "s",
      time_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 7: impute date to mid, time to last, preserve = FALSE ----
test_that("impute_dtc Test 7: impute date to mid, time to last, preserve = FALSE", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-15T00:00:00",
    "2019-06-30T00:00:00",
    "2019-06-30T00:00:00",
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    "2020-07-15T00:00:00"
  )
  expect_equal(
    imputes <- impute_dtc(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "MID",
      time_imputation = "FIRST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 8: impute date to MID day/month, time to 00:00:00 ----
test_that("impute_dtc Test 8: impute date to MID day/month, time to 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-15T00:00:00",
    "2019-06-30T00:00:00",
    "2019-06-30T00:00:00",
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    "2020-07-15T00:00:00"
  )
  expect_equal(
    impute_dtc(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "MID"
    ),
    expected_output
  )
})

## Test 9: impute date to MID and preserve works ----
test_that("impute_dtc Test 9: impute date to MID and preserve works", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-15T00:00:00",
    "2019-06-30T00:00:00",
    "2019-06-07T00:00:00",
    "2003-12-15T00:15:18",
    "2003-12-15T13:00:19",
    "2020-07-15T00:00:00"

  )

  actual_output <- impute_dtc(
    dtc = input,
    highest_imputation = "M",
    date_imputation = "MID",
    preserve = TRUE
  )

  expect_equal(
    actual_output,
    expected_output
  )
})

## Test 10: impute to 01-01 and preserve works ----
test_that("impute_dtc Test 10: impute to 01-01 and preserve works", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-01T00:00:00",
    "2019-01-01T00:00:00",
    "2019-01-07T00:00:00",
    "2003-12-15T00:15:18",
    "2003-12-15T13:00:19",
    "2020-07-01T00:00:00"
  )

  actual_output <- impute_dtc(
    dtc = input,
    highest_imputation = "M",
    date_imputation = "01-01",
    preserve = TRUE
  )

  expect_equal(
    actual_output,
    expected_output
  )
})

## Test 11: min_dates parameter works ----
test_that("impute_dtc Test 11: min_dates parameter works", {
  expect_equal(
    impute_dtc(c("2020-12", "2020-11", NA_character_),
      min_dates = list(
        c(
          ymd_hms("2020-12-06T12:12:12"),
          NA,
          NA
        ),
        c(
          ymd_hms("2020-11-11T11:11:11"),
          ymd_hms("2020-11-11T11:11:11"),
          ymd_hms("2020-11-11T11:11:11")
        )
      ),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    c("2020-12-06T12:12:12", "2020-11-11T11:11:11", "2020-11-11T11:11:11")
  )
})

## Test 12: max_dates parameter works ----
test_that("impute_dtc Test 12: max_dates parameter works", {
  expect_equal(
    impute_dtc(c("2020-12", "2020-11", NA_character_),
      max_dates = list(
        c(ymd_hms("2020-12-06T12:12:12"), NA, ymd_hms("2020-09-13T08:30:00")),
        c(ymd(""), ymd("2020-11-11"), ymd(""))
      ),
      highest_imputation = "Y",
      date_imputation = "last"
    ),
    c("2020-12-06T12:12:12", "2020-11-11T23:59:59", "2020-09-13T08:30:00")
  )
})
