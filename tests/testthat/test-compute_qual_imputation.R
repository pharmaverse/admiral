## Test 33: inform about deprecated arguments ----
test_that("compute_qual_imputation Test 33: inform about deprecated arguments", {
  withr::local_options(lifecycle_verbosity = "warning")
  expect_snapshot({
    compute_qual_imputation("<40", imputation_type = 1)
    compute_qual_imputation("<40", imputation_type = 2)
    compute_qual_imputation("<40", factor = 0.5)
  })
})

test_that("compute_qual_imputation Test 1:  Method 1 , > Qualifier Test", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = "strip"), 40.1)
})

test_that("compute_qual_imputation Test 2: Method 1, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = "strip"), 40.1)
})

test_that("compute_qual_imputation Test 3: Method 1, > Qualifier", {
  expect_equal(compute_qual_imputation(">40", imputation_type = "strip"), 40)
})

test_that("compute_qual_imputation Test 4: Method 1, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = "strip"), 40)
})

test_that("compute_qual_imputation Test 5: Method 1, < Qualifier", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = "strip"), 40.1)
})

test_that("compute_qual_imputation Test 6: Method 1, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = "strip"), 40.1)
})

test_that("compute_qual_imputation Test 7: Method 1, < Qualifier", {
  expect_equal(compute_qual_imputation("<40", imputation_type = "strip"), 40)
})

test_that("compute_qual_imputation Test 8: Method 1, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = "strip"), 40)
})

###############

test_that("compute_qual_imputation Test 9: Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = "strip", offset = 10), 50.1)
})

test_that("compute_qual_imputation Test 10: Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = "strip", offset = 10), 40.1)
})

test_that("compute_qual_imputation Test 11: Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">40", imputation_type = "strip", offset = 10), 50)
})

test_that("compute_qual_imputation Test 12: Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = "strip", offset = 10), 40)
})

test_that("compute_qual_imputation Test 13: Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = "strip", offset = 10), 30.1)
})

test_that("compute_qual_imputation Test 14: Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = "strip", offset = 10), 40.1)
})

test_that("compute_qual_imputation Test 15: Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<40", imputation_type = "strip", offset = 10), 30)
})

test_that("compute_qual_imputation Test 16: Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = "strip", offset = 10), 40)
})

###############

test_that("compute_qual_imputation Test 17: Method 2, > Qualifier", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = "shift"), 40.2)
})

test_that("compute_qual_imputation Test 18: Method 2, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = "shift"), 40.1)
})

test_that("compute_qual_imputation Test 19: Method 2, > Qualifier", {
  expect_equal(compute_qual_imputation(">40", imputation_type = "shift"), 41)
})

test_that("compute_qual_imputation Test 20: Method 2, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = "shift"), 40)
})

test_that("compute_qual_imputation Test 21: Method 2, < Qualifier", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = "shift"), 40.0)
})

test_that("compute_qual_imputation Test 22: Method 2, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = "shift"), 40.1)
})

test_that("compute_qual_imputation Test 23: Method 2, < Qualifier", {
  expect_equal(compute_qual_imputation("<40", imputation_type = "shift"), 39)
})

test_that("compute_qual_imputation Test 24: Method 2, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = "shift"), 40)
})

## Test 24a: documented vector examples ----
test_that("compute_qual_imputation Test 24a: documented vector examples", {
  values <- c("<40", ">40.1", "<=40", ">=40.1", "3")

  expect_equal(compute_qual_imputation(values), c(40, 40.1, 40, 40.1, 3))
  expect_equal(
    compute_qual_imputation(values, offset = 0.5),
    c(39.5, 40.6, 40, 40.1, 3)
  )
  expect_equal(
    compute_qual_imputation(values, imputation_type = "shift"),
    c(39, 40.2, 40, 40.1, 3)
  )
  expect_equal(
    compute_qual_imputation(values, imputation_type = "shift", offset = 100),
    c(39, 40.2, 40, 40.1, 3)
  )
  expect_equal(
    compute_qual_imputation(c(">40.1", ">40.10"), imputation_type = "shift"),
    c(40.2, 40.11)
  )
})

############

test_that("compute_qual_imputation Test 25: Method 2, No Qualifier", {
  expect_equal(compute_qual_imputation("40", imputation_type = "shift"), 40)
})

test_that("compute_qual_imputation Test 26: Method 2, No Qualifier, Character Value", {
  expect_equal(compute_qual_imputation("AB", imputation_type = "strip"), NA_real_)
})

test_that("compute_qual_imputation Test 27: Vector of Values", {
  expect_equal(compute_qual_imputation(c("3", ">30.2")), c(3.0, 30.2))
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

## Test 29: preserve valid results and decimal precision ----
test_that("compute_qual_imputation Test 29: preserve valid results and decimal precision", {
  values <- c("<40", ">40.1", "<=40", ">=40.1", "3", ">40.10", "<0.01", ">0", ">1.")
  expect_equal(compute_qual_imputation(values), c(40, 40.1, 40, 40.1, 3, 40.1, 0.01, 0, 1))
  expect_equal(
    compute_qual_imputation(values, "strip", offset = 0.5),
    c(39.5, 40.6, 40, 40.1, 3, 40.6, -0.49, 0.5, 1.5)
  )
  expect_equal(
    compute_qual_imputation(values, "shift"),
    c(39, 40.2, 40, 40.1, 3, 40.11, 0, 1, 2)
  )
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
    expect_error(compute_qual_imputation("<40", method), class = "compute_qual_imputation_method")
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

## Test 32: preserve deprecated calls ----
test_that("compute_qual_imputation Test 32: preserve deprecated calls", {
  withr::local_options(lifecycle_verbosity = "quiet")
  values <- c("<40", ">40.1", "<=40", ">=40.1", "3")
  expect_equal(compute_qual_imputation(values, 1, 0.5), c(39.5, 40.6, 40, 40.1, 3))
  expect_equal(compute_qual_imputation(values, 2, 100), c(39, 40.2, 40, 40.1, 3))
  expect_equal(compute_qual_imputation(values, factor = 0.5), c(39.5, 40.6, 40, 40.1, 3))
  expect_equal(compute_qual_imputation(values, 2L), c(39, 40.2, 40, 40.1, 3))
})
