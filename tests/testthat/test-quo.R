# quo_c ----
## Test 1: issues deprecation warning ----
test_that("quo_c Test 1: issues deprecation warning", {
  expect_warning(
    quo_c(quo(USUBJID), quo(STUDYID)),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 2: `quo_c` works in concatenating and indexing quosures ----
test_that("quo_c Test 2: `quo_c` works in concatenating and indexing quosures", {
  x <- quo(USUBJID)
  y <- quo(STUDYID)

  expect_equal(
    expected = quo(USUBJID),
    object = quo_c(x, NULL, y)[[1]]
  )
  expect_equal(
    expected = quo(STUDYID),
    object = quo_c(x, NULL, y)[[2]]
  )
})

## Test 3: `quo_c` returns error if non-quosures are input ----
test_that("quo_c Test 3: `quo_c` returns error if non-quosures are input", {
  USUBJID <- "01-701-1015" # nolint

  expect_error(
    object = quo_c(quo(USUBJID), USUBJID)
  )
})

# expr_c ----
## Test 4: concatenating and indexing expressions ----
test_that("expr_c Test 4: concatenating and indexing expressions", {
  x <- expr(USUBJID)
  y <- expr(STUDYID)

  expect_equal(
    expected = expr(USUBJID),
    object = expr_c(x, NULL, y)[[1]]
  )
  expect_equal(
    expected = expr(STUDYID),
    object = expr_c(x, NULL, y)[[2]]
  )
})

## Test 5: returns error if non-expressions are input ----
test_that("expr_c Test 5: returns error if non-expressions are input", {
  expect_error(
    object = expr_c(expr(USUBJID), mean)
  )
})

# quo_not_missing ----
## Test 6: issues deprecation warning ----
test_that("quo_not_missing Test 6: issues deprecation warning", {
  test_fun <- function(x) {
    x <- enquo(x)
    !isTRUE(quo_not_missing(x))
  }
  expect_warning(
    test_fun(my_variable),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 7: `quo_not_missing` throws an Error if missing argument ----
test_that("quo_not_missing Test 7: `quo_not_missing` throws an Error if missing argument", {
  test_fun <- function(x) {
    x <- enquo(x)
    isTrue(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})

# replace_values_by_names ----
## Test 8: names of quosures replace value ----
test_that("replace_values_by_names Test 8: names of quosures replace value", {
  z <- exprs(USUBJID, STUDYID)

  z_noname <- replace_values_by_names(z)

  names(z) <- c("Unique Subject Identifier", "Study Identifier")
  z_named <- replace_values_by_names(z)

  expect_equal(
    expected = expr(USUBJID),
    object = z_noname[[1]]
  )
  expect_equal(
    expected = expr(STUDYID),
    object = z_noname[[2]]
  )

  expect_equal(
    expected = expr(`Unique Subject Identifier`),
    object = z_named[[1]]
  )
  expect_equal(
    expected = expr(`Study Identifier`),
    object = z_named[[2]]
  )
})

## Test 9: warning if quosures argument is used ----
test_that("replace_values_by_names Test 9: warning if quosures argument is used", {
  expect_warning(
    replace_values_by_names(quosures = rlang::quos(STUDYID, USUBJID)),
    class = "lifecycle_warning_deprecated"
  )
})

# replace_symbol_in_quo ----
## Test 10: error if called ----
test_that("replace_symbol_in_quo Test 10: error if called", {
  expect_error(
    replace_symbol_in_quo(),
    class = "lifecycle_error_deprecated"
  )
})

# replace_symbol_in_expr ----
## Test 11: symbol is replaced ----
test_that("replace_symbol_in_expr Test 11: symbol is replaced", {
  expect_equal(
    expected = expr(AVAL.join),
    object = replace_symbol_in_expr(
      expr(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 12: partial match is not replaced ----
test_that("replace_symbol_in_expr Test 12: partial match is not replaced", {
  expect_equal(
    expected = expr(AVALC),
    object = replace_symbol_in_expr(
      expr(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 13: symbol in expression is replaced ----
test_that("replace_symbol_in_expr Test 13: symbol in expression is replaced", {
  expect_equal(
    expected = expr(desc(AVAL.join)),
    object = replace_symbol_in_expr(
      expr(desc(AVAL)),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

# add_suffix_to_vars ----
## Test 14: with single variable ----
test_that("add_suffix_to_vars Test 14: with single variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL),
      suffix = ".join"
    )
  )
})

## Test 15: with more than one variable ----
test_that("add_suffix_to_vars Test 15: with more than one variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
