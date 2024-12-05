test_that("my_first_fcn works correctly", {
  result <- my_first_fcn()
  expect_type(result, "character") # Test that the output is a string
  expect_equal(result, "Welcome to the admiral family!") # Test the exact content
})
