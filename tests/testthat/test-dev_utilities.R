# arg_name ----
## Test 1: arg_name works ----
test_that("arg_name Test 1: arg_name works", {
  expect_equal(arg_name(sym("a")), "a")
  expect_equal(arg_name(call("enquo", sym("a"))), "a")
  expect_error(arg_name("a"), "Could not extract argument name from")
})

# convert_dtm_to_dtc ----
## Test 2: works if dtm is in correct format ----
test_that("convert_dtm_to_dtc Test 2: works if dtm is in correct format", {
  expect_equal(
    convert_dtm_to_dtc(as.POSIXct("2022-04-05 15:34:07 UTC")),
    "2022-04-05T15:34:07"
  )
})

## Test 3: Error is thrown if dtm is not in correct format ----
test_that("convert_dtm_to_dtc Test 3: Error is thrown if dtm is not in correct format", {
  expect_error(
    convert_dtm_to_dtc("2022-04-05T15:26:14"),
    "lubridate::is.instant(dtm) is not TRUE",
    fixed = TRUE
  )
})

test_that("filter_if Test 1 : Input is returned as is if filter is NULL", {
  library(tibble)

  input <- tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expected_output <- input

  expect_dfs_equal(
    expected_output,
    filter_if(input, quo(NULL)),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("filter_if Test 2 : Input is filtered if filter is not NULL", {
  library(tibble)

  input <- tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expected_output <- input[1L, ]

  expect_dfs_equal(
    expected_output,
    filter_if(input, quo(VSTESTCD == "WEIGHT")),
    keys = c("USUBJID", "VSTESTCD")
  )
})
