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

# impute_dtc_dtm ----
## Test 1: default: no date imputation, time part set to 00:00:00 ----
test_that("impute_dtc_dtm Test 1: default: no date imputation, time part set to 00:00:00", {
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
  expect_equal(impute_dtc_dtm(dtc = input), expected_output)
})

## Test 2: no date imputation, min and sec imputed with 59 ----
test_that("impute_dtc_dtm Test 2: no date imputation, min and sec imputed with 59", {
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
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "23:59:59"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "LAST"
    ),
    expected_output
  )
})

## Test 3: impute date to first, time to 00:00:00 ----
test_that("impute_dtc_dtm Test 3: impute date to first, time to 00:00:00", {
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
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "01-01"
    ),
    expected_output
  )
})

## Test 4: impute day to last, time to 23:59:59 ----
test_that("impute_dtc_dtm Test 4: impute day to last, time to 23:59:59", {
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
    impute_dtc_dtm(
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
test_that("impute_dtc_dtm Test 5: impute date to last, time to 23:59:59 and preserve = TRUE", {
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
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "LAST",
      time_imputation = "LAST",
      preserve = TRUE
    ),
    expected_output
  )
})

## Test 6: no date imputation, impute second to 59 ----
test_that("impute_dtc_dtm Test 6: no date imputation, impute second to 59", {
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
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "s",
      time_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 7: impute date to mid, time to last ----
test_that("impute_dtc_dtm Test 7: impute date to mid, time to last", {
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
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "MID",
      time_imputation = "FIRST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 8: min_dates parameter works ----
test_that("impute_dtc_dtm Test 8: min_dates parameter works", {
  expect_equal(
    impute_dtc_dtm(c("2020-12", "2020-11", NA_character_),
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

## Test 9: max_dates parameter works ----
test_that("impute_dtc_dtm Test 9: max_dates parameter works", {
  expect_equal(
    impute_dtc_dtm(c("2020-12", "2020-11", NA_character_),
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

# impute_dtc_dt ----
input <- c(
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 10: default: no date imputation ----
test_that("impute_dtc_dt Test 10: default: no date imputation", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    NA_character_,
    NA_character_
  )
  expect_equal(impute_dtc_dt(dtc = input), expected_output)
})

## Test 11: impute date to first ----
test_that("impute_dtc_dt Test 11: impute date to first", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-01",
    "2019-01-01",
    "2019-01-01"
  )

  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "01-01"
    ),
    expected_output
  )
})

## Test 12: impute day to last ----
test_that("impute_dtc_dt Test 12: impute day to last", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-28",
    NA_character_,
    NA_character_
  )
  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "D",
      date_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 13: impute date to last and preserve = TRUE ----
test_that("impute_dtc_dt Test 13: impute date to last and preserve = TRUE", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-28",
    "2019-12-31",
    "2019-12-07"
  )
  expect_equal(
    imputes <- impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "LAST",
      preserve = TRUE
    ),
    expected_output
  )
})


## Test 14: impute date to mid ----
test_that("impute_dtc_dt Test 14: impute date to mid", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-15",
    "2019-06-30",
    "2019-06-30"
  )
  expect_equal(
    imputes <- impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "mid"
    ),
    expected_output
  )
})

## Test 15: min_dates parameter works ----
test_that("impute_dtc_dt Test 15: min_dates parameter works", {
  expect_equal(
    impute_dtc_dt(
      c("2020-12", "2020-11", NA_character_),
      min_dates = list(c(ymd("2020-12-06"),
                         NA,
                         NA),
                       c(
                         ymd("2020-11-11"),
                         ymd("2020-11-11"),
                         ymd("2020-11-11")
                       )),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    c("2020-12-06", "2020-11-11", "2020-11-11")
  )
})

## Test 16: max_dates parameter works ----
test_that("impute_dtc_dt Test 16: max_dates parameter works", {
  expect_equal(
    impute_dtc_dt(c("2020-12", "2020-11", NA_character_),
                   max_dates = list(
                     c(ymd("2020-12-06"), NA, ymd("2020-09-13")),
                     c(ymd(""), ymd("2020-11-11"), ymd(""))
                   ),
                   highest_imputation = "Y",
                   date_imputation = "last"
    ),
    c("2020-12-06", "2020-11-11", "2020-09-13")
  )
})

# convert_dtc_to_dtm ----
## Test 17: Convert a complete -- DTC into a date time object ----
test_that("convert_dtc_to_dtm Test 17: Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_dtc_to_dtm("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

# convert_dtc_to_dt ----
inputdtc <- c(
  "2019-07-18T15:25:52",
  "2019-07-18"
)

## Test 18: Convert a complete -- DTC into a date object ----
test_that("convert_dtc_to_dt Test 18: Convert a complete -- DTC into a date object", {
  expected_output <- c(
    as.Date("2019-07-18"),
    as.Date("2019-07-18")
  )
  expect_equal(
    convert_dtc_to_dt(dtc = inputdtc),
    expected_output
  )
})

# convert_date_to_dtm
## Test 19: Convert a complete -- DTC into a date time object ----
test_that("convert_date_to_dtm Test 19: Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

## Test 20: Impute incomplete -- DTC into a date time object ----
test_that("convert_date_to_dtm Test 20: Impute incomplete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18", time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

## Test 21: Convert -- DT into a date time object ----
test_that("convert_date_to_dtm Test 21: Convert -- DT into a date time object", {
  expect_equal(
    convert_date_to_dtm(as.Date("2019-07-18"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

##  Test 22: Keep -- DTM as the original date time object ----
test_that("convert_date_to_dtm Test 22: Keep -- DTM as the original date time object", {
  expect_equal(
    convert_date_to_dtm(ymd_hms("2019-07-18T15:25:52"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

# compute_dtf ----

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

## Test 23: compute DTF ----
test_that("compute_dtf Test 23: compute DTF", {
  expected_output <- c(
    NA_character_,
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

# compute_tmf ----
## Test 24: compute TMF ----
test_that("compute_tmf Test 24: compute TMF", {
  input_dtc <- c(
    "2019-07-18T15:25:52",
    "2019-07-18T15:25",
    "2019-07-18T15",
    "2019-07-18",
    "2019-02",
    "2019",
    "2019---07",
    "2003-12-15T-:15:18",
    "2003-12-15T13:-:19",
    "2020-07--T00:00",
    "2020-07--T00:00:00"
  )
  input_dtm <- c(
    as.POSIXct("2019-07-18T15:25:52"),
    as.POSIXct("2019-07-18T15:25:00"),
    as.POSIXct("2019-07-18T15:00:00"),
    as.POSIXct("2019-07-18"),
    as.POSIXct("2019-02-01"),
    as.POSIXct("2019-01-01"),
    as.POSIXct("2019-01-01"),
    as.POSIXct("2003-12-15T23:15:18"),
    as.POSIXct("2003-12-15T13:59:19"),
    as.POSIXct("2020-07-31T00:00:59"),
    as.POSIXct("2020-07-31T00:00:59")
  )
  expected_output <- c(
    NA_character_,
    "S",
    "M",
    "H",
    "H",
    "H",
    "H",
    "H",
    "M",
    "S",
    NA_character_
  )

  expect_equal(
    compute_tmf(
      dtc = input_dtc,
      dtm = input_dtm
    ),
    expected_output
  )
})
