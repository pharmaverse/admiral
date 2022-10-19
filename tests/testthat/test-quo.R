# quo_c ----
## Test 1: `quo_c` works in concatenating and indexing quosures ----
test_that("quo_c Test 1: `quo_c` works in concatenating and indexing quosures", {
  x <- quo(USUBJID)
  y <- quo(STUDYID)

  expect_equal(quo_c(x, NULL, y)[[1]], quo(USUBJID))
  expect_equal(quo_c(x, NULL, y)[[2]], quo(STUDYID))
})

## Test 2: `quo_c` returns error if non-quosures are input ----
test_that("quo_c Test 2: `quo_c` returns error if non-quosures are input", {
  USUBJID <- "STUDY-1001"

  expect_error(quo_c(quo(USUBJID), USUBJID))
})

# quo_not_missing ----
## Test 3: `quo_not_missing` returns TRUE if no missing argument ----
test_that("quo_not_missing Test 3: `quo_not_missing` returns TRUE if no missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_true(test_fun(my_variable))
})

## Test 4: `quo_not_missing` throws and Error if missing argument ----
test_that("quo_not_missing Test 4: `quo_not_missing` throws and Error if missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})

# replace_values_by_names ----
## Test 5: names of quosures replace value ----
test_that("replace_values_by_names Test 5: names of quosures replace value", {
  x <- quo(USUBJID)
  y <- quo(STUDYID)
  z <- quo_c(x, y)
  names(z) <- c("Unique Subject Identifier", "Study Identifier")
  z2 <- replace_values_by_names(z)

  expect_equal(z2[[1]], quo(`Unique Subject Identifier`))
  expect_equal(z2[[2]], quo(`Study Identifier`))
})

# replace_symbol_in_quo ----
## Test 6: symbol is replaced ----
test_that("replace_symbol_in_quo Test 6: symbol is replaced", {
  expect_equal(
    expected = quo(AVAL.join),
    object = replace_symbol_in_quo(
      quo(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 7: partial match is not replaced ----
test_that("replace_symbol_in_quo Test 7: partial match is not replaced", {
  expect_equal(
    expected = quo(AVALC),
    object = replace_symbol_in_quo(
      quo(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 8: symbol in expression is replaced ----
test_that("replace_symbol_in_quo Test 8: symbol in expression is replaced", {
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
## Test 9: with single variable ----
test_that("add_suffix_to_vars Test 9: with single variable", {
  expect_equal(
    expected = vars(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(vars(ADT, desc(AVAL), AVALC), vars = vars(AVAL), suffix = ".join")
  )
})

## Test 10: with more than one variable ----
test_that("add_suffix_to_vars Test 10: with more than one variable", {
  expect_equal(
    expected = vars(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      vars(ADT, desc(AVAL), AVALC),
      vars = vars(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
