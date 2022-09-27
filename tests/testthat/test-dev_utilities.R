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

# replace_symbol_in_quo ----
## Test 4: symbol is replaced ----
test_that("replace_symbol_in_quo Test 4: symbol is replaced", {
  expect_equal(
    expected = quo(AVAL.join),
    object = replace_symbol_in_quo(
      quo(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 5: partial match is not replaced ----
test_that("replace_symbol_in_quo Test 5: partial match is not replaced", {
  expect_equal(
    expected = quo(AVALC),
    object = replace_symbol_in_quo(
      quo(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 6: symbol in expression is replaced ----
test_that("replace_symbol_in_quo Test 6: symbol in expression is replaced", {
  expect_equal(
    expected = quo(desc(AVAL.join)),
    object = replace_symbol_in_quo(
      quo(desc(AVAL)),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

# add_suffix_to_vars ----
## Test 7: with single variable ----
test_that("add_suffix_to_vars Test 7: with single variable", {
  expect_equal(
    expected = vars(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(vars(ADT, desc(AVAL), AVALC), vars = vars(AVAL), suffix = ".join")
  )
})

## Test 8: with more than one variable ----
test_that("add_suffix_to_vars Test 8: with more than one variable", {
  expect_equal(
    expected = vars(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      vars(ADT, desc(AVAL), AVALC),
      vars = vars(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
