# my_first_fcn ----

## Test 1: input is Null ----
test_that("my_first_fcn Test 1: Null as input", {
  input <- NULL
  expected_output <- "Incorrect input"
  expect_equal(my_first_fcn(input), expected_output)
})

## Test 2: input is non-null ----
test_that("my_first_fcn Test 2: non-null input", {
  input <- "Y"
  expected_output <- "Welcome to the admiral family!"
  expect_equal(my_first_fcn(input), expected_output)
})
