test_that("compute_qual_imputation Test 1:  Method 1 , > Qualifier Test", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = 1), 40.1)
})

test_that("compute_qual_imputation Test 2: Method 1, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = 1), 40.1)
})

test_that("compute_qual_imputation Test 3: Method 1, > Qualifier", {
  expect_equal(compute_qual_imputation(">40", imputation_type = 1), 40)
})

test_that("compute_qual_imputation Test 4: Method 1, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = 1), 40)
})

test_that("compute_qual_imputation Test 5: Method 1, < Qualifier", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = 1), 40.1)
})

test_that("compute_qual_imputation Test 6: Method 1, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = 1), 40.1)
})

test_that("compute_qual_imputation Test 7: Method 1, < Qualifier", {
  expect_equal(compute_qual_imputation("<40", imputation_type = 1), 40)
})

test_that("compute_qual_imputation Test 8: Method 1, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = 1), 40)
})

###############

test_that("compute_qual_imputation Test 9: Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = 1, factor = 10), 50.1)
})

test_that("compute_qual_imputation Test 10: Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = 1, factor = 10), 40.1)
})

test_that("compute_qual_imputation Test 11: Method 1, > Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">40", imputation_type = 1, factor = 10), 50)
})

test_that("compute_qual_imputation Test 12: Method 1, >= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = 1, factor = 10), 40)
})

test_that("compute_qual_imputation Test 13: Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = 1, factor = 10), 30.1)
})

test_that("compute_qual_imputation Test 14: Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = 1, factor = 10), 40.1)
})

test_that("compute_qual_imputation Test 15: Method 1, < Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<40", imputation_type = 1, factor = 10), 30)
})

test_that("compute_qual_imputation Test 16: Method 1, <= Qualifier, Factor=10", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = 1, factor = 10), 40)
})

###############

test_that("compute_qual_imputation Test 17: Method 2, > Qualifier", {
  expect_equal(compute_qual_imputation(">40.1", imputation_type = 2), 40.2)
})

test_that("compute_qual_imputation Test 18: Method 2, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40.1", imputation_type = 2), 40.1)
})

test_that("compute_qual_imputation Test 19: Method 2, > Qualifier", {
  expect_equal(compute_qual_imputation(">40", imputation_type = 2), 41)
})

test_that("compute_qual_imputation Test 20: Method 2, >= Qualifier", {
  expect_equal(compute_qual_imputation(">=40", imputation_type = 2), 40)
})

test_that("compute_qual_imputation Test 21: Method 2, < Qualifier", {
  expect_equal(compute_qual_imputation("<40.1", imputation_type = 2), 40.0)
})

test_that("compute_qual_imputation Test 22: Method 2, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40.1", imputation_type = 2), 40.1)
})

test_that("compute_qual_imputation Test 23: Method 2, < Qualifier", {
  expect_equal(compute_qual_imputation("<40", imputation_type = 2), 39)
})

test_that("compute_qual_imputation Test 24: Method 2, <= Qualifier", {
  expect_equal(compute_qual_imputation("<=40", imputation_type = 2), 40)
})

############

test_that("compute_qual_imputation Test 25: Method 2, No Qualifier", {
  expect_equal(compute_qual_imputation("40", imputation_type = 2), 40)
})

test_that("compute_qual_imputation Test 26: Method 2, No Qualifier, Character Value", {
  expect_equal(compute_qual_imputation("AB", imputation_type = 1), NA_real_)
})
