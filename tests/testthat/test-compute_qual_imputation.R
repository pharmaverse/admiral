## Test 1: strip, > decimal qualifier ----
test_that("compute_qual_imputation Test 1: strip, > decimal qualifier", {
  expect_equal(compute_qual_imputation(c(">40.1", ">40.10", "3"), "strip"), c(40.1, 40.1, 3))
})

## Test 2: strip, >= decimal qualifier ----
test_that("compute_qual_imputation Test 2: strip, >= decimal qualifier", {
  expect_equal(compute_qual_imputation(c(">=40.1", ">=40.10", "3"), "strip"), c(40.1, 40.1, 3))
})

## Test 3: strip, > integer qualifier ----
test_that("compute_qual_imputation Test 3: strip, > integer qualifier", {
  expect_equal(compute_qual_imputation(c(">40", ">0", "3"), "strip"), c(40, 0, 3))
})

## Test 4: strip, >= integer qualifier ----
test_that("compute_qual_imputation Test 4: strip, >= integer qualifier", {
  expect_equal(compute_qual_imputation(c(">=40", ">=0", "3"), "strip"), c(40, 0, 3))
})

## Test 5: strip, < decimal qualifier ----
test_that("compute_qual_imputation Test 5: strip, < decimal qualifier", {
  expect_equal(compute_qual_imputation(c("<40.1", "<0.01", "3"), "strip"), c(40.1, 0.01, 3))
})

## Test 6: strip, <= decimal qualifier ----
test_that("compute_qual_imputation Test 6: strip, <= decimal qualifier", {
  expect_equal(compute_qual_imputation(c("<=40.1", "<=0.01", "3"), "strip"), c(40.1, 0.01, 3))
})

## Test 7: strip, < integer qualifier ----
test_that("compute_qual_imputation Test 7: strip, < integer qualifier", {
  expect_equal(compute_qual_imputation(c("<40", "<1", "3"), "strip"), c(40, 1, 3))
})

## Test 8: strip, <= integer qualifier ----
test_that("compute_qual_imputation Test 8: strip, <= integer qualifier", {
  expect_equal(compute_qual_imputation(c("<=40", "<=1", "3"), "strip"), c(40, 1, 3))
})

## Test 9: strip, > decimal qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 9: strip, > decimal qualifier, offset = 10", {
  expect_equal(
    compute_qual_imputation(c(">40.1", ">40.10", "3"), "strip", offset = 10),
    c(50.1, 50.1, 3)
  )
  expect_equal(compute_qual_imputation(c(">40.1", "<40", "3"), offset = 0.5), c(40.6, 39.5, 3))
})

## Test 10: strip, >= decimal qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 10: strip, >= decimal qualifier, offset = 10", {
  expect_equal(
    compute_qual_imputation(c(">=40.1", ">=40.10", "3"), "strip", offset = 10),
    c(40.1, 40.1, 3)
  )
})

## Test 11: strip, > integer qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 11: strip, > integer qualifier, offset = 10", {
  expect_equal(compute_qual_imputation(c(">40", ">0", "3"), "strip", offset = 10), c(50, 10, 3))
})

## Test 12: strip, >= integer qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 12: strip, >= integer qualifier, offset = 10", {
  expect_equal(compute_qual_imputation(c(">=40", ">=0", "3"), "strip", offset = 10), c(40, 0, 3))
})

## Test 13: strip, < decimal qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 13: strip, < decimal qualifier, offset = 10", {
  expect_equal(
    compute_qual_imputation(c("<40.1", "<0.01", "3"), "strip", offset = 10),
    c(30.1, -9.99, 3)
  )
})

## Test 14: strip, <= decimal qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 14: strip, <= decimal qualifier, offset = 10", {
  expect_equal(
    compute_qual_imputation(c("<=40.1", "<=0.01", "3"), "strip", offset = 10),
    c(40.1, 0.01, 3)
  )
})

## Test 15: strip, < integer qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 15: strip, < integer qualifier, offset = 10", {
  expect_equal(compute_qual_imputation(c("<40", "<1", "3"), "strip", offset = 10), c(30, -9, 3))
})

## Test 16: strip, <= integer qualifier, offset = 10 ----
test_that("compute_qual_imputation Test 16: strip, <= integer qualifier, offset = 10", {
  expect_equal(compute_qual_imputation(c("<=40", "<=1", "3"), "strip", offset = 10), c(40, 1, 3))
})

## Test 17: shift, > decimal qualifier and trailing zeros ----
test_that("compute_qual_imputation Test 17: shift, > decimal qualifier and trailing zeros", {
  expect_equal(compute_qual_imputation(c(">40.1", ">40.10", ">1."), "shift"), c(40.2, 40.11, 2))
})

## Test 18: shift, >= decimal qualifier ----
test_that("compute_qual_imputation Test 18: shift, >= decimal qualifier", {
  expect_equal(compute_qual_imputation(c(">=40.1", ">=40.10", "3"), "shift"), c(40.1, 40.1, 3))
})

## Test 19: shift, > integer qualifier ----
test_that("compute_qual_imputation Test 19: shift, > integer qualifier", {
  expect_equal(compute_qual_imputation(c(">40", ">0", "3"), "shift"), c(41, 1, 3))
})

## Test 20: shift, >= integer qualifier ----
test_that("compute_qual_imputation Test 20: shift, >= integer qualifier", {
  expect_equal(compute_qual_imputation(c(">=40", ">=0", "3"), "shift"), c(40, 0, 3))
})

## Test 21: shift, < decimal qualifier ----
test_that("compute_qual_imputation Test 21: shift, < decimal qualifier", {
  expect_equal(compute_qual_imputation(c("<40.1", "<0.01", "3"), "shift"), c(40, 0, 3))
})

## Test 22: shift, <= decimal qualifier ----
test_that("compute_qual_imputation Test 22: shift, <= decimal qualifier", {
  expect_equal(compute_qual_imputation(c("<=40.1", "<=0.01", "3"), "shift"), c(40.1, 0.01, 3))
})

## Test 23: shift, < integer qualifier ----
test_that("compute_qual_imputation Test 23: shift, < integer qualifier", {
  expect_equal(compute_qual_imputation(c("<40", "<1", "3"), "shift"), c(39, 0, 3))
})

## Test 24: shift, <= integer qualifier ----
test_that("compute_qual_imputation Test 24: shift, <= integer qualifier", {
  expect_equal(compute_qual_imputation(c("<=40", "<=1", "3"), "shift"), c(40, 1, 3))
})

## Test 25: shift, no qualifier ----
test_that("compute_qual_imputation Test 25: shift, no qualifier", {
  expect_equal(compute_qual_imputation(c("40", "-3", ".5"), "shift"), c(40, -3, 0.5))
})

## Test 26: strip, invalid character value ----
test_that("compute_qual_imputation Test 26: strip, invalid character value", {
  expect_equal(compute_qual_imputation(c("AB", NA_character_), "strip"), c(NA_real_, NA_real_))
})

## Test 27: default method preserves vector order ----
test_that("compute_qual_imputation Test 27: default method preserves vector order", {
  expect_equal(compute_qual_imputation(c("3", ">30.2", "<=40")), c(3, 30.2, 40))
})

## Test 28: reject malformed qualifiers silently ----
test_that("compute_qual_imputation Test 28: reject malformed qualifiers silently", {
  values <- c(
    "0>1", "-<0.9>", "< 40", " <40", "<40 ", "<40\n", "<<40", "><40",
    "=<40", "==40", "40<", "<", "<=", "<-40", ">.5", "1.2.3",
    "", " ", "AB", "1e3", "Inf", "NaN", NA_character_
  )
  for (method in c("strip", "shift")) {
    expect_no_warning(actual <- compute_qual_imputation(values, method))
    expect_identical(actual, rep(NA_real_, length(values)))
  }
  expect_equal(
    compute_qual_imputation(c("0>1", "-<0.9>", ">100.00", "<1"), "shift"),
    c(NA_real_, NA_real_, 100.01, 0)
  )
})

## Test 29: preserve unqualified, empty and missing results ----
test_that("compute_qual_imputation Test 29: preserve unqualified, empty and missing results", {
  expect_equal(compute_qual_imputation(c("-3", "+3", ".5", "-.5", "3.")), c(-3, 3, 0.5, -0.5, 3))
  expect_identical(compute_qual_imputation(character()), numeric())
  expect_identical(compute_qual_imputation(character(), "shift"), numeric())
  expect_identical(
    compute_qual_imputation(c(NA_character_, NA_character_), "shift"),
    c(NA_real_, NA_real_)
  )
})

## Test 30: validate argument types and method ----
test_that("compute_qual_imputation Test 30: validate argument types and method", {
  for (value in list(1, TRUE, factor("<40"), list("<40"), NULL)) {
    expect_error(compute_qual_imputation(value), class = "assert_character_vector")
  }
  for (method in list(
    "other", "1", NA_character_, character(), c("strip", "shift"), TRUE, NULL, 1 + 0i
  )) {
    expect_error(compute_qual_imputation("<40", method), class = "assert_character_scalar")
  }
  for (method in list(0, 3, 1.5, NA_real_, numeric(), c(1, 2), Inf)) {
    expect_error(
      compute_qual_imputation("<40", method),
      regexp = '"strip" or "shift"',
      class = "compute_qual_imputation_method"
    )
  }
  expect_snapshot(
    compute_qual_imputation("<40", imputation_type = 0),
    error = TRUE, cnd_class = TRUE
  )
})

## Test 31: validate offsets ----
test_that("compute_qual_imputation Test 31: validate offsets", {
  for (offset in list("1", TRUE, NULL, numeric(), c(1, 2), 1 + 0i)) {
    expect_error(compute_qual_imputation("<40", offset = offset), class = "assert_numeric_vector")
  }
  for (offset in list(NA_real_, NaN, Inf, -Inf, -1)) {
    expect_error(
      compute_qual_imputation("<40", offset = offset),
      class = "compute_qual_imputation_offset"
    )
  }
  expect_snapshot(
    compute_qual_imputation("<40", offset = -1),
    error = TRUE, cnd_class = TRUE
  )
  expect_error(compute_qual_imputation("<40", factor = "1"), class = "assert_numeric_vector")
  expect_snapshot(
    compute_qual_imputation("<40", factor = 1, offset = 1),
    error = TRUE, cnd_class = TRUE
  )
  expect_equal(compute_qual_imputation("<40", "shift", offset = 100), 39)
})

## Test 32: preserve deprecated calls and messages ----
test_that("compute_qual_imputation Test 32: preserve deprecated calls and messages", {
  withr::local_options(lifecycle_verbosity = "warning")
  values <- c("<40", ">40.1", "<=40", ">=40.1", "3")
  expect_snapshot({
    expect_equal(compute_qual_imputation(values, imputation_type = 1), c(40, 40.1, 40, 40.1, 3))
    expect_equal(compute_qual_imputation(values, imputation_type = 2), c(39, 40.2, 40, 40.1, 3))
    expect_equal(compute_qual_imputation(values, 1, 0.5), c(39.5, 40.6, 40, 40.1, 3))
    expect_equal(compute_qual_imputation(values, 2, 100), c(39, 40.2, 40, 40.1, 3))
    expect_equal(compute_qual_imputation(values, factor = 0.5), c(39.5, 40.6, 40, 40.1, 3))
    expect_equal(compute_qual_imputation(values, 2L), c(39, 40.2, 40, 40.1, 3))
  })
})
