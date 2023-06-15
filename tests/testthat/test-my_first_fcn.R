## Test 1: hello_admiral() function greets ----
test_that("my_first_fcn Test 1: hello_admiral function greets", {
  expect_message(
    hello_admiral(),
    "Welcome to the admiral family!"
  )
})
