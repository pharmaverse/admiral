# my_first_fcn ----
## Test 1: returns NULL ----
test_that("my_first_fcn Test 1: returns NULL", {
  expect_null(my_first_fcn())
})

## Test 2: produces a message ----
test_that("my_first_fcn Test 2: produces a message", {
  expect_message(my_first_fcn(), "^Welcome to the admiral family!\\n$")
})
