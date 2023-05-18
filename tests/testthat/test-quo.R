# quo_c ----
## Test 1: issues deprecation warning ----
test_that("quo_c Test 1: issues deprecation warning", {
  expect_error(
    quo_c(quo(USUBJID), quo(STUDYID)),
    class = "lifecycle_error_deprecated"
  )
})

## Test 2: `quo_c` works in concatenating and indexing quosures ----
test_that("quo_c Test 2: `quo_c` works in concatenating and indexing quosures", {
  x <- quo(USUBJID)
  y <- quo(STUDYID)

  expect_error(
    quo_c(x, NULL, y)[[1]],
    class = "lifecycle_error_deprecated"
  )
  expect_error(
    object = quo_c(x, NULL, y)[[2]],
    class = "lifecycle_error_deprecated"
  )
})

## Test 3: `quo_c` returns error if non-quosures are input ----
test_that("quo_c Test 3: `quo_c` returns error if non-quosures are input", {
  USUBJID <- "01-701-1015" # nolint

  expect_error(
    object = quo_c(quo(USUBJID), USUBJID),
    class = "lifecycle_error_deprecated"
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

## Test 5: concatenating named list of expressions ----
test_that("expr_c Test 5: concatenating named list of expressions", {
  expect_equal(
    expected = exprs(PARAMCD = "DOSE", PARAMN = 1),
    object = expr_c(exprs(PARAMCD = "DOSE", PARAMN = 1, NULL))
  )
})

## Test 6: concatenating list and single expression ----
test_that("expr_c Test 6: concatenating list and single expression", {
  expect_equal(
    expected = exprs(PARAMCD, PARAM, AVAL),
    object = expr_c(exprs(PARAMCD, PARAM), expr(AVAL))
  )
})

## Test 7: returns error if non-expressions are input ----
test_that("expr_c Test 7: returns error if non-expressions are input", {
  expect_error(
    object = expr_c(expr(USUBJID), mean)
  )
})

# quo_not_missing ----
## Test 8: issues deprecation warning ----
test_that("quo_not_missing Test 8: issues deprecation warning", {
  test_fun <- function(x) {
    x <- enquo(x)
    !isTRUE(quo_not_missing(x))
  }
  expect_error(
    test_fun(my_variable),
    class = "lifecycle_error_deprecated"
  )
})

## Test 9: `quo_not_missing` throws an Error if missing argument ----
test_that("quo_not_missing Test 9: `quo_not_missing` throws an Error if missing argument", {
  test_fun <- function(x) {
    x <- enquo(x)
    isTrue(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})

# replace_values_by_names ----
## Test 10: names of expressions replace value ----
test_that("replace_values_by_names Test 10: names of expressions replace value", {
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

## Test 11: names of argument is NULL ----
test_that("replace_values_by_names Test 11: names of argument is NULL", {
  z <- exprs(USUBJID, STUDYID)
  names(z) <- NULL

  expect_equal(
    expected = z,
    object = replace_values_by_names(z)
  )
})

## Test 12: warning if quosures argument is used ----
test_that("replace_values_by_names Test 12: warning if quosures argument is used", {
  expect_error(
    replace_values_by_names(quosures = rlang::quos(STUDYID, USUBJID)),
    class = "lifecycle_error_deprecated"
  )
})

# replace_symbol_in_quo ----
## Test 13: error if called ----
test_that("replace_symbol_in_quo Test 13: error if called", {
  expect_error(
    replace_symbol_in_quo(),
    class = "lifecycle_error_deprecated"
  )
})

# replace_symbol_in_expr ----
## Test 14: symbol is replaced ----
test_that("replace_symbol_in_expr Test 14: symbol is replaced", {
  expect_equal(
    expected = expr(AVAL.join),
    object = replace_symbol_in_expr(
      expr(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 15: partial match is not replaced ----
test_that("replace_symbol_in_expr Test 15: partial match is not replaced", {
  expect_equal(
    expected = expr(AVALC),
    object = replace_symbol_in_expr(
      expr(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 16: symbol in expression is replaced ----
test_that("replace_symbol_in_expr Test 16: symbol in expression is replaced", {
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
## Test 17: with single variable ----
test_that("add_suffix_to_vars Test 17: with single variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL),
      suffix = ".join"
    )
  )
})

## Test 18: with more than one variable ----
test_that("add_suffix_to_vars Test 18: with more than one variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
