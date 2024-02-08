# my_first_fcn
## Test 1: Welcome message while not the first time ----
test_that(
  "my_first_fcn Test 1: Welcome message while not the first time",
  {
    expect_message(my_first_fcn(first_time = FALSE), "Welcome again to the admiral family!")
  }
)

## Test 2: Welcome message while this is the first time ----
test_that(
  "my_first_fcn Test 2: Welcome message while this is the first time",
  {
    expect_message(my_first_fcn(first_time = TRUE), "Welcome to the admiral family!")
  }
)
