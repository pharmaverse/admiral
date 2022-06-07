test_that("compute_qualifier_imputation Test 1:  Method 1 , > Qualifier Test", {
  expect_equal(compute_qualifier_imputation(">40.1", imputation_type = 1), 40.1)
})

test_that("compute_qualifier_imputation Test 2: Method 1, >= Qualifier", {
  expect_equal(compute_qualifier_imputation(">=40.1", imputation_type = 1), 40.1)
})

test_that("compute_qualifier_imputation Test 3: Method 1, > Qualifier", {
  expect_equal(compute_qualifier_imputation(">40", imputation_type = 1), 40)
})

test_that("compute_qualifier_imputation Test 4: Method 1, >= Qualifier", {
  expect_equal(compute_qualifier_imputation(">=40", imputation_type = 1), 40)
})

test_that("compute_qualifier_imputation Test 5: Method 1, < Qualifier", {
  expect_equal(compute_qualifier_imputation("<40.1", imputation_type = 1), 40.1)
})

test_that("Method 1, <= Qualifier", {
  expect_equal(compute_qualifier_imputation("<=40.1", imputation_type = 1), 40.1)
})

test_that("Method 1, < Qualifier", {
  expect_equal(compute_qualifier_imputation("<40", imputation_type = 1), 40)
})

test_that("Method 1, <= Qualifier", {
  expect_equal(compute_qualifier_imputation("<=40", imputation_type = 1), 40)
})

###############

test_that("Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation(">40.1", imputation_type = 1, factor = 10), 50.1)
})

test_that("Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation(">=40.1", imputation_type = 1, factor = 10), 40.1)
})

test_that("Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation(">40", imputation_type = 1, factor = 10), 50)
})

test_that("Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation(">=40", imputation_type = 1, factor = 10), 40)
})

test_that("Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation("<40.1", imputation_type = 1, factor = 10), 30.1)
})

test_that("Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation("<=40.1", imputation_type = 1, factor = 10), 40.1)
})

test_that("Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation("<40", imputation_type = 1, factor = 10), 30)
})

test_that("Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qualifier_imputation("<=40", imputation_type = 1, factor = 10), 40)
})

###############

test_that("Method 2, > Qualifier", {
  expect_equal(compute_qualifier_imputation(">40.1", imputation_type = 2), 40.2)
})

test_that("Method 2, >= Qualifier", {
  expect_equal(compute_qualifier_imputation(">=40.1", imputation_type = 2), 40.1)
})

test_that("Method 2, > Qualifier", {
  expect_equal(compute_qualifier_imputation(">40", imputation_type = 2), 41)
})

test_that("Method 2, >= Qualifier", {
  expect_equal(compute_qualifier_imputation(">=40", imputation_type = 2), 40)
})

test_that("Method 2, < Qualifier", {
  expect_equal(compute_qualifier_imputation("<40.1", imputation_type = 2), 40.0)
})

test_that("Method 2, <= Qualifier", {
  expect_equal(compute_qualifier_imputation("<=40.1", imputation_type = 2), 40.1)
})

test_that("Method 2, < Qualifier", {
  expect_equal(compute_qualifier_imputation("<40", imputation_type = 2), 39)
})

test_that("Method 2, <= Qualifier", {
  expect_equal(compute_qualifier_imputation("<=40", imputation_type = 2), 40)
})

############

test_that("Method 2, No Qualifier", {
  expect_equal(compute_qualifier_imputation("40", imputation_type = 2), 40)
})

test_that("Method 2, No Qualifier, Character Value", {
  expect_equal(compute_qualifier_imputation("AB", imputation_type = 1), NA_real_)
})
