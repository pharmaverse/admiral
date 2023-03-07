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

# filter_if ----
## Test 4: Input is returned as is if filter is NULL ----
test_that("filter_if Test 4: Input is returned as is if filter is NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01",    "WEIGHT",       80.9,
    "P01",    "HEIGHT",      189.2
  )

  expected_output <- input

  expect_dfs_equal(
    expected_output,
    filter_if(input, expr(NULL)),
    keys = c("USUBJID", "VSTESTCD")
  )
})

## Test 5: Input is filtered if filter is not NULL ----
test_that("filter_if Test 5: Input is filtered if filter is not NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01",    "WEIGHT",       80.9,
    "P01",    "HEIGHT",      189.2
  )

  expected_output <- input[1L, ]

  expect_dfs_equal(
    expected_output,
    filter_if(input, expr(VSTESTCD == "WEIGHT")),
    keys = c("USUBJID", "VSTESTCD")
  )
})

# contains_vars ----
## Test 6: returns TRUE for valid arguments ----
test_that("contains_vars Test 6: returns TRUE for valid arguments", {
  expect_true(contains_vars(exprs(USUBJID, PARAMCD)))
})

## Test 7: returns TRUE for valid arguments ----
test_that("contains_vars Test 7: returns TRUE for valid arguments", {
  expect_error(contains_vars(USUBJID))
})
# vars2chr ----
## Test 8: returns character vector ----
test_that("vars2chr Test 8: returns character vector", {
  expected <- c("STUDYID", "USUBJID")
  names(expected) <- c("", "")
  expect_equal(vars2chr(exprs(STUDYID, USUBJID)), expected)
})

## Test 9: warning if quosures argument is used ----
test_that("vars2chr Test 9: warning if quosures argument is used", {
  expect_warning(
    vars2chr(quosures = rlang::quos(STUDYID, USUBJID)),
    class = "lifecycle_warning_deprecated"
  )
})
