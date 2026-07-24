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

# get_dt_dtm_range ----
## Test 13: get_dt_dtm_range correctly imputes date ranges ----
test_that("get_dt_dtm_range Test 13: get_dt_dtm_range correctly imputes date ranges", {

  # Validate date range imputation for complete dates and partial dates
  expect_equal(
    get_dt_dtm_range(
      c("2020-02-29", "2021-03", "--08-15", "2022", "2022--", "---29"),
      create_datetime = FALSE
    ),
    list(
      lower = c(
        "2020-02-29", "2021-03-01", "0000-01-01", "2022-01-01",
        "2022-01-01", "0000-01-01"
      ),
      upper = c(
        "2020-02-29", "2021-03-31", "9999-12-31", "2022-12-31",
        "2022-12-31", "9999-12-31"
      )
    )
  )

  # Check leap year February
  expect_equal(
    get_dt_dtm_range(c("2020-02", "2021-02"), create_datetime = FALSE),
    list(
      lower = c("2020-02-01", "2021-02-01"),
      upper = c("2020-02-29", "2021-02-28")
    )
  )

  # Validate datetime range imputation
  expect_equal(
    get_dt_dtm_range(
      c(
        "2020-02-29T12:00", "2021-03T14:30", "2021-03--T14",
        "--03-15T14:00", "2022--T12:34"
      ),
      create_datetime = TRUE
    ),
    list(
      lower = c(
        "2020-02-29T12:00:00", "2021-03-01T00:00:00", "2021-03-01T00:00:00",
        "0000-01-01T00:00:00", "2022-01-01T00:00:00"
      ),
      upper = c(
        "2020-02-29T12:00:59", "2021-03-31T23:59:59", "2021-03-31T23:59:59",
        "9999-12-31T23:59:59", "2022-12-31T23:59:59"
      )
    )
  )

  # Edge case: empty input
  expect_equal(get_dt_dtm_range(character(0), create_datetime = FALSE), character(0))

  # Validate imputation with invalid date formats (warnings)
  invalid_dtc <- c("invalid-date", "2021-13-40")
  expect_snapshot(
    get_dt_dtm_range(invalid_dtc, create_datetime = FALSE)
  )
})

## Test 14: lower_bounds ----
test_that("get_dt_dtm_range Test 14: lower_bounds", {
  expect_equal(
    get_dt_dtm_range(
      c("2020-02-20", "2021-03", "2021-03", NA),
      lower_bounds = list(
        c(ymd("2020-02-10"), ymd("2021-03-03"), ymd("2021-04-01"), ymd("2022-01-01")),
        c(NA, ymd("2021-03-05"), NA, ymd("2021-10-11"))
      ),
      create_datetime = FALSE
    ),
    list(
      lower = c("2020-02-20", "2021-03-05", "2021-03-01", "2022-01-01"),
      upper = c("2020-02-20", "2021-03-31", "2021-03-31", "9999-12-31")
    )
  )
})

## Test 15: upper_bounds ----
test_that("get_dt_dtm_range Test 15: upper_bounds", {
  expect_equal(
    get_dt_dtm_range(
      c("2020-02-09", "2021-03", "2021-03", NA),
      upper_bounds = list(
        c(ymd("2020-02-05"), ymd("2021-03-23"), ymd("2021-04-01"), ymd("2022-01-21")),
        c(NA, ymd("2021-03-15"), NA, ymd("2021-10-11"))
      ),
      create_datetime = FALSE
    ),
    list(
      lower = c("2020-02-09", "2021-03-01", "2021-03-01", "0000-01-01"),
      upper = c("2020-02-09", "2021-03-15", "2021-03-31", "2021-10-11")
    )
  )
})


# get_highest_imputation_level ----
## Test 16: correctly determine highest level ----
test_that("get_highest_imputation_level Test 16: correctly determine highest level", {
  expect_equal(get_highest_imputation_level("Y", FALSE), dt_level("Y"))
  expect_equal(get_highest_imputation_level("Y", TRUE), dtm_level("Y"))
})

# get_imputation_targets ----
## Test 17: correctly generate imputation targets ----
test_that("get_imputation_targets Test 17: correctly generate imputation targets", {
  # Date tests
  partial_date <- list(year = "2020", month = "03", day = NA_character_)

  # When date_imputation = "first"
  expect_equal(
    get_imputation_targets(partial_date, date_imputation = "first"),
    list(year = "0000", month = "01", day = "01")
  )

  # When date_imputation = "mid"
  expect_equal(
    get_imputation_targets(partial_date, date_imputation = "mid"),
    list(year = "xxxx", month = "06", day = "15")
  )

  partial_with_na_month <- list(year = "2020", month = NA_character_, day = NA_character_)
  expect_equal(
    get_imputation_targets(partial_with_na_month, date_imputation = "mid"),
    list(year = "xxxx", month = "06", day = "30")
  )

  # When date_imputation = "last"
  expect_equal(
    get_imputation_targets(partial_date, date_imputation = "last"),
    list(year = "9999", month = "12", day = "28")
  )

  # When date_imputation = "mm-dd"
  expect_equal(
    get_imputation_targets(partial_date, date_imputation = "06-15"),
    list(year = "xxxx", month = "06", day = "15")
  )

  # Datetime tests
  partial_datetime <- list(
    year = "2020", month = "03", day = NA_character_,
    hour = "12", minute = NA_character_, second = NA_character_
  )

  # When date_imputation = "first" and time_imputation = "first"
  expect_equal(
    get_imputation_targets(partial_datetime, date_imputation = "first", time_imputation = "first"),
    list(year = "0000", month = "01", day = "01", hour = "00", minute = "00", second = "00")
  )

  # When date_imputation = "last" and time_imputation = "last"
  expect_equal(
    get_imputation_targets(partial_datetime, date_imputation = "last", time_imputation = "last"),
    list(year = "9999", month = "12", day = "28", hour = "23", minute = "59", second = "59")
  )

  # When date_imputation = "mid" and time_imputation = "first"
  expect_equal(
    get_imputation_targets(partial_datetime, date_imputation = "mid", time_imputation = "first"),
    list(year = "xxxx", month = "06", day = "15", hour = "00", minute = "00", second = "00")
  )

  # When date_imputation = "06-15" and time_imputation = "12:34:56"
  expect_equal(
    get_imputation_targets(partial_datetime,
      date_imputation = "06-15",
      time_imputation = "12:34:56"
    ),
    list(year = "xxxx", month = "06", day = "15", hour = "12", minute = "34", second = "56")
  )

  # Test error when time_imputation is NULL for datetime
  expect_error(
    get_imputation_targets(partial_datetime, date_imputation = "first", time_imputation = NULL),
    "As `partial` is datetime, `time_imputation` is expected."
  )
})


# adjust_last_day_imputation ----
## Test 18: correctly adjust last day imputation ----
test_that("adjust_last_day_imputation Test 18: correctly adjust last day imputation", {
  expect_equal(
    adjust_last_day_imputation(
      imputed_dtc = "2021-03-01",
      partial = list(year = "2021", month = "03", day = NA_character_)
    ),
    "2021-03-31"
  )
  expect_equal(
    adjust_last_day_imputation(
      imputed_dtc = "2021-03-01T00:00:00",
      partial = list(
        year = "2021", month = "03", day = NA_character_,
        hour = NA_character_, minute = NA_character_, second = NA_character_
      )
    ),
    "2021-03-31T00:00:00"
  )
})

# impute_date_time ----
## Test 19: correctly imputes missing values ----
test_that("impute_date_time Test 19: correctly imputes missing values", {
  partial <- list(year = "2020", month = NA_character_, day = NA_character_)
  target <- list(year = "2020", month = "01", day = "01")
  expect_equal(
    impute_date_time(partial, target),
    list(year = "2020", month = "01", day = "01")
  )
})

## Test 20: gives error if partial and target differ. ----
test_that("impute_date_time Test 20: gives error if partial and target differ.", {
  expect_error(
    admiral:::impute_date_time(
      partial = list(year = "2020", month = "05", day = NA_character_),
      target = list(year = "2020", day = "05", hour = "12")
    ),
    regexp = "Names of `partial` and `target` do not match."
  )
})

# format_imputed_dtc ----
## Test 21: correctly format imputed date time ----
test_that("format_imputed_dtc Test 21: correctly format imputed date time", {
  imputed <- list(
    year = "2020", month = "01", day = "01",
    hour = "12", minute = "00", second = "00"
  )
  expect_equal(format_imputed_dtc(imputed), "2020-01-01T12:00:00")
  expect_equal(format_imputed_dtc(imputed[1:3]), "2020-01-01")
})

# propagate_na_values ----
## Test 22: correctly propagate NA values in DateTime ----
test_that("propagate_na_values Test 22: correctly propagate NA values in DateTime", {
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

# get_partialdatetime ----
## Test 23: get_partialdatetime parses complete datetime ----
test_that("get_partialdatetime Test 23: get_partialdatetime parses complete datetime", {
  dtc <- "2020-12-31T23:59:59"
  expected <- list(
    year = "2020",
    month = "12",
    day = "31",
    hour = "23",
    minute = "59",
    second = "59"
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 24: get_partialdatetime parses partial datetime ----
test_that("get_partialdatetime Test 24: get_partialdatetime parses partial datetime", {
  dtc <- "2020-12-31T23:59"
  expected <- list(
    year = "2020",
    month = "12",
    day = "31",
    hour = "23",
    minute = "59",
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 25: get_partialdatetime parses date-only input ----
test_that("get_partialdatetime Test 25: get_partialdatetime parses date-only input", {
  dtc <- "2020-12-31"
  expected <- list(
    year = "2020",
    month = "12",
    day = "31",
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 26: get_partialdatetime parses year-only input ----
test_that("get_partialdatetime Test 26: get_partialdatetime parses year-only input", {
  dtc <- "2020"
  expected <- list(
    year = "2020",
    month = NA_character_,
    day = NA_character_,
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 27: create_datetime = FALSE only returns date parts ----
test_that("get_partialdatetime Test 27: create_datetime = FALSE only returns date parts", {
  dtc <- "2020-03-15T12:34:56"
  expected <- list(
    year = "2020",
    month = "03",
    day = "15"
  )
  result <- admiral:::get_partialdatetime(dtc, FALSE)
  expect_identical(result, expected)
})

## Test 28: handle partial date with missing components ----
test_that("get_partialdatetime Test 28: handle partial date with missing components", {
  dtc <- "2020--"
  expected <- list(
    year = "2020",
    month = NA_character_,
    day = NA_character_,
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 29: handle partial datetime with some time components ----
test_that("get_partialdatetime Test 29: handle partial datetime with some time components", {
  dtc <- "2020-03-15T12"
  expected <- list(
    year = "2020",
    month = "03",
    day = "15",
    hour = "12",
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 30: handle completely missing components ----
test_that("get_partialdatetime Test 30: handle completely missing components", {
  dtc <- "---T::"
  expected <- list(
    year = NA_character_,
    month = NA_character_,
    day = NA_character_,
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 31: handle empty string input ----
test_that("get_partialdatetime Test 31: handle empty string input", {
  dtc <- ""
  expected <- list(
    year = NA_character_,
    month = NA_character_,
    day = NA_character_,
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

## Test 32: handle NA input ----
test_that("get_partialdatetime Test 32: handle NA input", {
  dtc <- NA_character_
  expected <- list(
    year = NA_character_,
    month = NA_character_,
    day = NA_character_,
    hour = NA_character_,
    minute = NA_character_,
    second = NA_character_
  )
  result <- admiral:::get_partialdatetime(dtc, TRUE)
  expect_identical(result, expected)
})

# is_partial_datetime ----
## Test 33: correctly identifies datetime and date partials ----
test_that("is_partial_datetime Test 33: correctly identifies datetime and date partials", {
  # Test with a full datetime
  partial_datetime <- list(
    year = "2023", month = "05", day = "15",
    hour = "14", minute = "30", second = "00"
  )
  expect_true(is_partial_datetime(partial_datetime))

  # Test with a date only, partially NA
  partial_date <- list(year = "2023", month = NA_character_, day = "15")
  expect_false(is_partial_datetime(partial_date))

  # Test with an invalid partial (missing components)
  partial_invalid <- list(year = "2023", month = "05", hour = "14")
  expect_snapshot(is_partial_datetime(partial_invalid),
    error = TRUE
  )

  # Test with an empty list
  expect_snapshot(is_partial_datetime(list()),
    error = TRUE
  )

  # Test with extra components
  partial_extra <- list(
    year = "2023", month = "05", day = "15",
    hour = "14", minute = "30", second = "00",
    millisecond = "500"
  )
  expect_true(is_partial_datetime(partial_extra))

  # Test with only time components
  partial_time_only <- list(hour = "14", minute = "30", second = "00")
  expect_snapshot(is_partial_datetime(partial_time_only),
    error = TRUE
  )
})

# assert_time_imputation ----
## Test 34: gives error when input not a valid format ----
test_that("assert_time_imputation Test 34: gives error when input not a valid format", {
  expect_error(
    admiral:::assert_time_imputation(c("25:00:00"), "H"),
    regexp = paste0(
      '`time_imputation` must be one of "first", "last" or time specified as',
      ' "hh:mm:ss": e.g. "12:00:00"'
    )
  )
})
