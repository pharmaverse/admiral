# my_first_fcn ----
## Test 1: <Explanation of the test> ----
test_that("multiplication works", {
  expected_output <- "Welcome to the admiral family!"

  actual_output <- my_first_fcn()

  expect_equal(expected_output, actual_output)
})
