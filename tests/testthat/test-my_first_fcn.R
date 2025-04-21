# my_first_fcn ----
## Test 1: Using the default `exclamation` argument returns 'Welcome to the admiral family!' ----
test_that("my_first_fcn Test 1: Default argument returns 'Welcome to the admiral family!'", {

  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family!\\n"
    )

})

## Test 2: Setting argument `exclamation` to `TRUE` returns 'Welcome to the admiral family!' ----
test_that("my_first_fcn Test 2: Argument `exclamation=TRUE` returns 'Welcome to the admiral family!'", {

  expect_message(
    hello_admiral(exclamation = TRUE),
    "^Welcome to the admiral family!\\n"
  )

})

## Test 3: Setting argument `exclamation` to `FALSE` returns 'Welcome to the admiral family!' ----
test_that("my_first_fcn Test 3: Argument `exclamation=FALSE` returns 'Welcome to the admiral family!'", {

  expect_message(
    hello_admiral(exclamation = FALSE),
    "^Welcome to the admiral family.\\n"
  )

})
