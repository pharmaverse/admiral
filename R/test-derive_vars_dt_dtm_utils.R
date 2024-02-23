# convert_date_to_dtm ----
## Test 1: Convert a complete -- DTC into a date time object ----
test_that("convert_date_to_dtm  Test 1: Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:52")
  )
})

## Test 2: Impute incomplete -- DTC into a date time object ----
test_that("convert_date_to_dtm  Test 2: Impute incomplete -- DTC into a date time object", {
  expect_equal(
    convert_date_to_dtm("2019-07-18", time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

## Test 3: Convert -- DT into a date time object ----
test_that("convert_date_to_dtm  Test 3: Convert -- DT into a date time object", {
  expect_equal(
    convert_date_to_dtm(as.Date("2019-07-18"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T23:59:59")
  )
})

## Test4 4: Keep -- DTM as the original date time object ----
test_that("convert_date_to_dtm  Test 4: Keep -- DTM as the original date time object", {
  expect_equal(
    convert_date_to_dtm(ymd_hms("2019-07-18T15:25:52"), time_imputation = "23:59:59"),
    ymd_hms("2019-07-18T15:25:52")
  )
})
