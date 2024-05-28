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
