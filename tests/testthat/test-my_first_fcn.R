# anything param types ----

## Test 1: NULL works as parameter ----
test_that("my_first_fcn Test 1: NULL works as parameter", {
  input <- NULL

  expected_output <- "Welcome to the admiral family!"

  expect_equal(welcome_fun(input), expected_output)
})



## Test 2: double works as parameter ----
test_that("my_first_fcn Test 2: double works as parameter", {
  input <- 3

  expected_output <- "Welcome to the admiral family!"

  expect_equal(welcome_fun(input), expected_output)
})


## Test 3: dataframe works as parameter ----
test_that("my_first_fcn Test 3: dataframe works as parameter", {
  input <- tibble::tribble(
    ~x,  ~y,
    "a", 1:3,
    "b", 4:6
  )

  expected_output <- "Welcome to the admiral family!"

  expect_equal(welcome_fun(input), expected_output)
})
