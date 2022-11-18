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

## Test 3: impute month and day to first, time to 00:00:00 ----
test_that("impute_dtc_dtm Test 3: impute month and day to first, time to 00:00:00", {
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
      date_imputation = "last",
      time_imputation = "last"
    ),
    expected_output
  )
})

## Test 5: impute month, day to last, time to 23:59:59, preserve = TRUE ----
test_that("impute_dtc_dtm Test 5: impute month, day to last, time to 23:59:59, preserve = TRUE", {
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
      date_imputation = "last",
      time_imputation = "last",
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

## Test 7: impute month and day to mid, time to first ----
test_that("impute_dtc_dtm Test 7: impute month and day to mid, time to first", {
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
      date_imputation = "mid",
      time_imputation = "first",
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
    impute_dtc_dtm(c("2020-12", "2020-11", NA_character_, "2020-02-02"),
      max_dates = list(
        c(ymd_hms("2020-12-06T12:12:12"), NA, ymd_hms("2020-09-13T08:30:00"), NA),
        c(ymd(""), ymd("2020-11-11"), ymd(""), ymd("2020-02-02"))
      ),
      highest_imputation = "Y",
      date_imputation = "last",
      time_imputation = "last"
    ),
    c("2020-12-06T12:12:12", "2020-11-11T23:59:59", "2020-09-13T08:30:00", "2020-02-02T23:59:59")
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

## Test 11: impute month and day to first ----
test_that("impute_dtc_dt Test 11: impute month and day to first", {
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

## Test 13: impute month and day to last and preserve = TRUE ----
test_that("impute_dtc_dt Test 13: impute month and day to last and preserve = TRUE", {
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


## Test 14: impute month and day to mid ----
test_that("impute_dtc_dt Test 14: impute month and day to mid", {
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
      min_dates = list(
        c(
          ymd("2020-12-06"),
          NA,
          NA
        ),
        c(
          ymd("2020-11-11"),
          ymd("2020-11-11"),
          ymd("2020-11-11")
        )
      ),
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
    as.POSIXct(NA_character_),
    as.POSIXct(NA_character_),
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
    NA_character_,
    NA_character_,
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

## Test 25: throws ERROR when ignore_seconds_flag  = T and seconds are present ----
test_that("compute_tmf Test 25: throws ERROR when ignore_seconds_flag  = T and seconds are present", { # nolint
  expect_error(
    compute_tmf(
      dtc = c("2020-11-11T11:11:11", "2020-11-11T11:11"),
      dtm = ymd_hms(c(
        "2020-11-11T11:11:11", "2020-11-11T11:11:00"
      )),
      ignore_seconds_flag = TRUE
    ),
    regexp = "Seconds detected in data while ignore_seconds_flag is invoked"
  )
})

## Test 26: ignore_seconds_flag  = TRUE ----
test_that("compute_tmf Test 26: ignore_seconds_flag  = TRUE", {
  expect_equal(
    compute_tmf(
      dtc = c("2020-11-11T11:11", "2020-11-11T11"),
      dtm = ymd_hms(c(
        "2020-11-11T11:11:00", "2020-11-11T11:00:00"
      )),
      ignore_seconds_flag = TRUE
    ),
    c(NA_character_, "M")
  )
})

# derive_vars_dt ----

date <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 27: default behavior ----
test_that("derive_vars_dt Test 27: default behavior", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18",          as.Date("2019-07-18"),
    "2019-02",             as.Date(NA),
    "2019",                as.Date(NA),
    "2019---07",           as.Date(NA)
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

## Test 28: no date imputation, add DTF ----
test_that("derive_vars_dt Test 28: no date imputation, add DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,                ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18",          as.Date("2019-07-18"), NA_character_,
    "2019-02",             as.Date(NA),           NA_character_,
    "2019",                as.Date(NA),           NA_character_,
    "2019---07",           as.Date(NA),           NA_character_
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    flag_imputation = "date",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

## Test 29: date imputed to first, auto DTF ----
test_that("derive_vars_dt Test 29: date imputed to first, auto DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,                ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18",          as.Date("2019-07-18"), NA_character_,
    "2019-02",             as.Date("2019-02-01"), "D",
    "2019",                as.Date("2019-01-01"), "M",
    "2019---07",           as.Date("2019-01-01"), "M"
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 30: date imputed to last, no DTF ----
test_that("derive_vars_dt Test 30: date imputed to last, no DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18",          as.Date("2019-07-18"),
    "2019-02",             as.Date("2019-02-28"),
    "2019",                as.Date("2019-12-31"),
    "2019---07",           as.Date("2019-12-31")
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    flag_imputation = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

# derive_vars_dtm ----

input <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 31: default behavior ----
test_that("derive_vars_dtm Test 31: default behavior", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), "H",
    "2019-02",             ymd_hms(NA),                    NA_character_,
    "2019",                ymd_hms(NA),                    NA_character_,
    "2019---07",           ymd_hms(NA),                    NA_character_
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 32: date imputed to first, auto DTM/TMF ----
test_that("derive_vars_dtm Test 32: date imputed to first, auto DTF/TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTDTF,       ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02",             ymd_hms("2019-02-01T00:00:00"), "D",           "H",
    "2019",                ymd_hms("2019-01-01T00:00:00"), "M",           "H",
    "2019---07",           ymd_hms("2019-01-01T00:00:00"), "M",           "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 33: date and time imputed to last, no DTF/TMF ----
test_that("derive_vars_dtm Test 33: date and time imputed to last, no DTF/TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDTM,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"),
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:59"),
    "2019-07-18T15",       ymd_hms("2019-07-18T15:59:59"),
    "2019-07-18",          ymd_hms("2019-07-18T23:59:59"),
    "2019-02",             ymd_hms("2019-02-28T23:59:59"),
    "2019",                ymd_hms("2019-12-31T23:59:59"),
    "2019---07",           ymd_hms("2019-12-31T23:59:59")
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "LAST",
    time_imputation = "LAST",
    flag_imputation = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 34: date and time imputed to last, DTF only ----
test_that("derive_vars_dtm Test 34: date and time imputed to last, DTF only", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDTM,                        ~AENDTF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:59"), NA_character_,
    "2019-07-18T15",       ymd_hms("2019-07-18T15:59:59"), NA_character_,
    "2019-07-18",          ymd_hms("2019-07-18T23:59:59"), NA_character_,
    "2019-02",             ymd_hms("2019-02-28T23:59:59"), "D",
    "2019",                ymd_hms("2019-12-31T23:59:59"), "M",
    "2019---07",           ymd_hms("2019-12-31T23:59:59"), "M"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    flag_imputation = "date"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 35: date imputed to MID, time to first, TMF only ----
test_that("derive_vars_dtm Test 35: date imputed to MID, time to first, TMF only", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), "H",
    "2019-02",             ymd_hms("2019-02-15T00:00:00"), "H",
    "2019",                ymd_hms("2019-06-30T00:00:00"), "H",
    "2019---07",           ymd_hms("2019-06-07T00:00:00"), "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "mid",
    flag_imputation = "time",
    preserve = TRUE
  )

  expect_dfs_equal(
    base = expected_output,
    comp = actual_output,
    keys = c("XXSTDTC")
  )
})

## Test 36: No re-derivation is done if --DTF variable already exists ----
test_that("derive_vars_dtm Test 36: No re-derivation is done if --DTF variable already exists", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTDTF,       ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02",             ymd_hms("2019-02-01T00:00:00"), "D",           "H",
    "2019",                ymd_hms("2019-01-01T00:00:00"), "MD",          "H",
    "2019---07",           ymd_hms("2019-01-01T00:00:00"), "M",           "H"
  ) %>%
    select(XXSTDTC, ASTDTF, everything())

  actual_output <- expect_message(
    derive_vars_dtm(
      mutate(input, ASTDTF = c(NA, NA, NA, NA, "D", "MD", "M")),
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    regexp = "^The .* variable is already present in the input dataset and will not be re-derived."
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 37: max_dates parameter works as expected ----
test_that("derive_vars_dtm Test 37: max_dates parameter works as expected", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,    ~ASTDTM,                        ~ASTDTF, ~ASTTMF,
    "2019-02",   ymd_hms("2019-02-10T00:00:00"), "D",     "H",
    "2019",      ymd_hms("2019-02-10T00:00:00"), "M",     "H",
    "2019---07", ymd_hms("2019-02-10T00:00:00"), "M",     "H"
  ) %>%
    mutate(DCUTDT = ymd_hms("2019-02-10T00:00:00"))

  actual_output <- derive_vars_dtm(
    select(expected_output, XXSTDTC, DCUTDT),
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    max_dates = vars(DCUTDT)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("XXSTDTC")
  )
})

input_secs <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)
