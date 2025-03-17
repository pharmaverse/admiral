test_that("my_first_fcn returns correct output", {
  # Test 1: Check return value
  expect_equal(my_first_fcn(), "Welcome to the admiral family!")

  # Test 2: Verify message output
  expect_message(my_first_fcn(), "Welcome to the admiral family!")
})
