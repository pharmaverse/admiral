# Test for the my_first_fcn function
test_that("my_first_fcn displays the correct welcome message", {
  # Test that the function displays the correct message
  expect_message(
    my_first_fcn(),
    "Welcome to the admiral family!"
  )
})
