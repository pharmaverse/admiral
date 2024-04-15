# my_new_func ----
## Test 1: This is a test for the onboarding dummy issue 1839
## to print "Welcome to the admiral family!" ----

test_that("my_new_func Test 1: Prints without input", {
  expect_output(
    my_first_fcn(),
    "Welcome to the admiral family!"
  )
})

test_that("my_new_func Test 1: Prints regardless of input", {
  expect_output(
    my_first_fcn(x),
    "Welcome to the admiral family!"
  )
})
