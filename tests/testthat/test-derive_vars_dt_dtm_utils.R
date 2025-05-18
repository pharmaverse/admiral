# convert_date_to_dtm ----
## Test 1: Convert a complete -- DTC into a date time object ----
test_that("convert_date_to_dtm Test 1: Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

## Test 2: Impute incomplete -- DTC into a date time object ----
test_that("convert_date_to_dtm Test 2: Impute incomplete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18", time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

## Test 3: Convert -- DT into a date time object ----
test_that("convert_date_to_dtm Test 3: Convert -- DT into a date time object", {
  expect_equal(
    convert_date_to_dtm(as.Date("2019-07-18"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

## Test 4: Keep -- DTM as the original date time object ----
test_that("convert_date_to_dtm Test 4: Keep -- DTM as the original date time object", {
  expect_equal(
    convert_date_to_dtm(ymd_hms("2019-07-18T15:25:52"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

# get_imputation_target_date ----
## Test 5: get correct target for missing dates ----
test_that("get_imputation_target_date Test 5: get correct target for missing dates", {
  expect_equal(
    get_imputation_target_date("first", NA),
    list(year = "0000", month = "01", day = "01")
  )
})


## Test 6: get correct target for missing dates ----
test_that("get_imputation_target_date Test 6: get correct target for missing dates", {
  expect_equal(
    get_imputation_target_date("mid", "04"),
    list(year = "xxxx", month = "06", day = "15")
  )
})

## Test 7: get correct target for missing dates ----
test_that("get_imputation_target_date Test 7: get correct target for missing dates", {
  expect_equal(
    get_imputation_target_date("mid", NA),
    list(year = "xxxx", month = "06", day = "30")
  )
})


## Test 8: get correct target for missing dates ----
test_that("get_imputation_target_date Test 8: get correct target for missing dates", {
  expect_equal(
    get_imputation_target_date("last", NA),
    list(year = "9999", month = "12", day = "28")
  )
})

# get_imputation_target_time ----
## Test 9: get correct target for missing times (first) ----
test_that("get_imputation_target_time Test 9: get correct target for missing times (first)", {
  expect_equal(
    get_imputation_target_time("first"),
    list(hour = "00", minute = "00", second = "00")
  )
})

## Test 10: get correct target for missing times (last) ----
test_that("get_imputation_target_time Test 10: get correct target for missing times (last)", {
  expect_equal(
    get_imputation_target_time("last"),
    list(hour = "23", minute = "59", second = "59")
  )
})

## Test 11: get correct target for missing times (time) ----
test_that("get_imputation_target_time Test 11: get correct target for missing times (time)", {
  expect_equal(
    get_imputation_target_time("22:04:35"),
    list(hour = "22", minute = "04", second = "35")
  )
})

## Test 12: get correct target for missing times (junk) ----
test_that("get_imputation_target_time Test 12: get correct target for missing times (junk)", {
  expect_equal(
    get_imputation_target_time("xxAyyAzz"),
    list(hour = "xx", minute = "yy", second = "zz")
  )
})

# Test for get_dt_dtm_range
## Test 13: get_dt_dtm_range correctly imputes date ranges ----
test_that("get_dt_dtm_range correctly imputes date ranges", {
  dtc <- c("2020-02-29", "2021-03", "--08-15", "2022", "2022--", "---29")
  dtc_datetimes <- c(
    "2020-02-29T12:00", "2021-03T14:30", "2021-03--T14",
    "--03-15T14:00", "2022--T12:34"
  )

  # Validate 'first' date imputation for complete dates and partial dates
  expect_equal(
    get_dt_dtm_range(dtc, date_imputation = "first"),
    c(
      "2020-02-29", "2021-03-01", "0000-08-15", "2022-01-01",
      "2022-01-01", "0000-01-29"
    )
  )

  # Validate 'last' date imputation for complete dates and partial dates
  expect_equal(
    get_dt_dtm_range(dtc, date_imputation = "last"),
    c(
      "2020-02-29", "2021-03-31", "9999-08-15", "2022-12-31",
      "2022-12-31", "9999-12-29"
    )
  )

  # Check leap year February
  expect_equal(
    get_dt_dtm_range(c("2020-02", "2021-02"), date_imputation = "last"),
    c("2020-02-29", "2021-02-28")
  )

  # Validate 'first' datetime imputation
  expect_equal(
    get_dt_dtm_range(dtc_datetimes,
      date_imputation = "first",
      time_imputation = "first"
    ),
    c(
      "2020-02-29T12:00:00", "2021-03-01T00:00:00", "2021-03-01T00:00:00",
      "0000-01-01T00:00:00", "2022-01-01T00:00:00"
    )
  )

  # Validate 'last' datetime imputation for partial datetime
  expect_equal(
    get_dt_dtm_range(dtc_datetimes,
      date_imputation = "first",
      time_imputation = "last"
    ),
    c(
      "2020-02-29T12:00:59", "2021-03-01T23:59:59", "2021-03-01T23:59:59",
      "0000-01-01T23:59:59", "2022-01-01T23:59:59"
    )
  )

  # Edge case: empty input
  expect_equal(get_dt_dtm_range(character(0), date_imputation = "first"), character(0))

  # Validate imputation with invalid date formats (warnings)
  invalid_dtc <- c("invalid-date", "2021-13-40")
  expect_warning(get_dt_dtm_range(invalid_dtc, date_imputation = "first"))
  expect_warning(get_dt_dtm_range(invalid_dtc, date_imputation = "last"))

  # Validate correct imputation with single year
  year_only <- c("2020", "2021")
  expect_equal(
    get_dt_dtm_range(year_only, date_imputation = "first"),
    c("2020-01-01", "2021-01-01")
  )
  expect_equal(
    get_dt_dtm_range(year_only, date_imputation = "last"),
    c("2020-12-31", "2021-12-31")
  )

  # Validate correct imputation with year-month only
  year_month <- c("2020-02", "2021-05")
  expect_equal(
    get_dt_dtm_range(year_month, date_imputation = "first"),
    c("2020-02-01", "2021-05-01")
  )
  expect_equal(
    get_dt_dtm_range(year_month, date_imputation = "last"),
    c("2020-02-29", "2021-05-31")
  )

  # Validate input with 'first' date and custom time imputation ('12:34:56')
  expect_equal(
    get_dt_dtm_range(dtc_datetimes,
      date_imputation = "first",
      time_imputation = "12:34:56"
    ),
    c(
      "2020-02-29T12:00:56", "2021-03-01T12:34:56", "2021-03-01T12:34:56",
      "0000-01-01T12:34:56", "2022-01-01T12:34:56"
    )
  )

  # Validate correct handling of date imputation 'first' with missing day
  dtc_partial <- c("2020-02", "2021-03T12")
  expect_equal(
    get_dt_dtm_range(dtc_partial,
      date_imputation = "first",
      time_imputation = "first"
    ),
    c("2020-02-01T00:00:00", "2021-03-01T00:00:00")
  )

  # Validate correct handling of date imputation 'last' with missing day
  expect_equal(
    get_dt_dtm_range(dtc_partial,
      date_imputation = "last",
      time_imputation = "last"
    ),
    c("2020-02-29T23:59:59", "2021-03-31T23:59:59")
  )
})

# Test for get_highest_imputation_level
## Test 14: correctly determine highest level ----
test_that("get_imputation_target_time Test 14: correctly determine highest level", {
  expect_equal(get_highest_imputation_level(FALSE, "Y"), dt_level("Y"))
  expect_equal(get_highest_imputation_level(TRUE, "Y"), dtm_level("Y"))
})

# Test for get_imputation_targets
## Test 15: correctly generate imputation targets ----
test_that("get_imputation_target_time Test 15: correctly generate imputation targets", {
  # When date_imputation = "first"
  partial <- list(year = "2020", month = "03", day = NA_character_)
  expect_equal(
    get_imputation_targets(partial,
      date_imputation = "first",
      time_imputation = "first",
      is_datetime = FALSE
    ),
    list(year = "0000", month = "01", day = "01")
  )
  expect_equal(
    get_imputation_targets(partial,
      date_imputation = "first",
      time_imputation = "first",
      is_datetime = TRUE
    ),
    c(
      list(year = "0000", month = "01", day = "01"),
      list(hour = "00", minute = "00", second = "00")
    )
  )

  # When date_imputation = "mid"
  expect_equal(
    get_imputation_targets(partial,
      date_imputation = "mid",
      time_imputation = "mid",
      is_datetime = FALSE
    ),
    list(year = "xxxx", month = "06", day = "15")
  )
  partial_with_na_month <- list(year = "2020", month = NA_character_, day = NA_character_)
  expect_equal(
    get_imputation_targets(partial_with_na_month,
      date_imputation = "mid",
      time_imputation = "mid",
      is_datetime = FALSE
    ),
    list(year = "xxxx", month = "06", day = "30")
  )

  # When date_imputation = "last"
  expect_equal(
    get_imputation_targets(partial,
      date_imputation = "last",
      time_imputation = "last",
      is_datetime = FALSE
    ),
    list(year = "9999", month = "12", day = "28")
  )

  # When date_imputation = "mm-dd"
  expect_equal(
    get_imputation_targets(partial,
      date_imputation = "06-15",
      time_imputation = NULL,
      is_datetime = FALSE
    ),
    list(year = "xxxx", month = "06", day = "15")
  )
})


# Test for adjust_last_day_imputation
## Test 16: correctly adjust last day imputation ----
test_that("get_imputation_target_time Test 16: correctly adjust last day imputation", {
  expect_equal(
    adjust_last_day_imputation(
      imputed_dtc = "2021-03-01",
      partial = list(day = NA_character_),
      is_datetime = FALSE
    ),
    "2021-03-31"
  )
  expect_equal(
    adjust_last_day_imputation(
      imputed_dtc = "2021-03-01T00:00:00",
      partial = list(day = NA_character_),
      is_datetime = TRUE
    ),
    "2021-03-31T00:00:00"
  )
})

# Test for impute_values
## Test 17: impute_values correctly imputes missing values ----
test_that("get_imputation_target_time Test 17: impute_values correctly imputes missing values", {
  partial <- list(year = "2020", month = NA_character_, day = NA_character_)
  target <- list(year = "2020", month = "01", day = "01")
  components <- c("year", "month", "day")
  expect_equal(
    impute_values(partial, target, components),
    list(year = "2020", month = "01", day = "01")
  )
})

# Test for format_imputed_dtc
## Test 18: correctly format imputed date time ----
test_that("get_imputation_target_time Test 18: correctly format imputed date time", {
  imputed <- list(
    year = "2020", month = "01", day = "01",
    hour = "12", minute = "00", second = "00"
  )
  expect_equal(format_imputed_dtc(imputed, TRUE), "2020-01-01T12:00:00")
  expect_equal(format_imputed_dtc(imputed[1:3], FALSE), "2020-01-01")
})

# Test for propagate_na_values
## Test 19: correctly propagate NA values in DateTime ----
test_that("get_imputation_target_time Test 19: correctly propagate NA values in DateTime", {
  partial <- list(
    year = "2020", month = NA_character_, day = "01",
    hour = "12", minute = NA_character_, second = "20"
  )
  expect_equal(
    propagate_na_values(partial),
    list(
      year = "2020", month = NA_character_, day = NA_character_,
      hour = NA_character_, minute = NA_character_, second = NA_character_
    )
  )
})

# Test for parse_partial_date_time
## Test 20: correctly parse partial date or datetime ----
test_that("get_imputation_target_time Test 20: correctly parse partial date or datetime", {
  dtc <- c("2020-02", "2021", "2020---10", "")
  expect_equal(
    parse_partial_date_time(dtc, FALSE),
    list(
      year = c("2020", "2021", "2020", NA_character_),
      month = c("02", NA_character_, NA_character_, NA_character_),
      day = c(NA_character_, NA_character_, "10", NA_character_)
    )
  )
  expect_equal(
    parse_partial_date_time(c(dtc, "2020-02-03T12"), TRUE),
    list(
      year = c("2020", "2021", "2020", NA_character_, "2020"),
      month = c("02", NA_character_, NA_character_, NA_character_, "02"),
      day = c(NA_character_, NA_character_, "10", NA_character_, "03"),
      hour = c(NA_character_, NA_character_, NA_character_, NA_character_, "12"),
      minute = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
      second = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    )
  )
})
