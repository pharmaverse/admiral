## Test 1: Verifying the message function works fine and prints out the message ----
test_that("Message printing works", {
  expected_output <- "Welcome to the admiral family"
  actual_output <- my_first_fcn()
  expect_equal(expected_output, actual_output)
})
