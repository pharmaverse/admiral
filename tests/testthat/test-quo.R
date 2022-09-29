# quo_not_missing ----
## Test 1: `quo_not_missing` returns TRUE if no missing argument ----
test_that("quo_not_missing Test 1: `quo_not_missing` returns TRUE if no missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_true(test_fun(my_variable))
})

## Test 2: `quo_not_missing` throws and Error if missing argument ----
test_that("quo_not_missing Test 2: `quo_not_missing` throws and Error if missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})

# replace_symbol_in_quo ----
## Test 3: symbol is replaced ----
test_that("replace_symbol_in_quo Test 3: symbol is replaced", {
  expect_equal(
    expected = quo(AVAL.join),
    object = replace_symbol_in_quo(
      quo(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 4: partial match is not replaced ----
test_that("replace_symbol_in_quo Test 4: partial match is not replaced", {
  expect_equal(
    expected = quo(AVALC),
    object = replace_symbol_in_quo(
      quo(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 5: symbol in expression is replaced ----
test_that("replace_symbol_in_quo Test 5: symbol in expression is replaced", {
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
## Test 6: with single variable ----
test_that("add_suffix_to_vars Test 6: with single variable", {
  expect_equal_tbl(
    expected = vars(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(vars(ADT, desc(AVAL), AVALC), vars = vars(AVAL), suffix = ".join")
  )
})

## Test 7: with more than one variable ----
test_that("add_suffix_to_vars Test 7: with more than one variable", {
  expect_equal_tbl(
    expected = vars(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      vars(ADT, desc(AVAL), AVALC),
      vars = vars(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
