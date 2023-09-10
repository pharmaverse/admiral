# expr_c ----
## Test 1: concatenating and indexing expressions ----
test_that("expr_c Test 1: concatenating and indexing expressions", {
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

## Test 2: concatenating named list of expressions ----
test_that("expr_c Test 2: concatenating named list of expressions", {
  expect_equal(
    expected = exprs(PARAMCD = "DOSE", PARAMN = 1),
    object = expr_c(exprs(PARAMCD = "DOSE", PARAMN = 1, NULL))
  )
})

## Test 3: concatenating list and single expression ----
test_that("expr_c Test 3: concatenating list and single expression", {
  expect_equal(
    expected = exprs(PARAMCD, PARAM, AVAL),
    object = expr_c(exprs(PARAMCD, PARAM), expr(AVAL))
  )
})

## Test 4: returns error if non-expressions are input ----
test_that("expr_c Test 4: returns error if non-expressions are input", {
  expect_error(
    object = expr_c(expr(USUBJID), mean)
  )
})

# replace_values_by_names ----
## Test 5: names of expressions replace value ----
test_that("replace_values_by_names Test 5: names of expressions replace value", {
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

## Test 6: names of argument is NULL ----
test_that("replace_values_by_names Test 6: names of argument is NULL", {
  z <- exprs(USUBJID, STUDYID)
  names(z) <- NULL

  expect_equal(
    expected = z,
    object = replace_values_by_names(z)
  )
})

## Test 7: warning if quosures argument is used ----
test_that("replace_values_by_names Test 7: warning if quosures argument is used", {
  expect_error(
    replace_values_by_names(quosures = rlang::quos(STUDYID, USUBJID)),
    class = "lifecycle_error_deprecated"
  )
})

# replace_symbol_in_expr ----
## Test 8: symbol is replaced ----
test_that("replace_symbol_in_expr Test 8: symbol is replaced", {
  expect_equal(
    expected = expr(AVAL.join),
    object = replace_symbol_in_expr(
      expr(AVAL),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 9: partial match is not replaced ----
test_that("replace_symbol_in_expr Test 9: partial match is not replaced", {
  expect_equal(
    expected = expr(AVALC),
    object = replace_symbol_in_expr(
      expr(AVALC),
      target = AVAL,
      replace = AVAL.join
    )
  )
})

## Test 10: symbol in expression is replaced ----
test_that("replace_symbol_in_expr Test 10: symbol in expression is replaced", {
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
## Test 11: with single variable ----
test_that("add_suffix_to_vars Test 11: with single variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL),
      suffix = ".join"
    )
  )
})

## Test 12: with more than one variable ----
test_that("add_suffix_to_vars Test 12: with more than one variable", {
  expect_equal(
    expected = exprs(ADT, desc(AVAL.join), AVALC.join),
    object = add_suffix_to_vars(
      exprs(ADT, desc(AVAL), AVALC),
      vars = exprs(AVAL, AVALC),
      suffix = ".join"
    )
  )
})
