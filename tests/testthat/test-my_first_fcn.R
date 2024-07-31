## Test my_first_fcn ---
test_that("my_first_fcn print correct output", {
  expect_message(my_first_fcn(), "Welcome to the admiral family!")
})
