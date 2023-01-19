# is_order_vars ----
## Test 1: returns error if input were created incorrectly ----
test_that("is_order_vars Test 1: returns error if input were created incorrectly", {
  expect_error(is_order_vars(STUDYID))
})

## Test 2: returns TRUE if input were created correctly ----
test_that("is_order_vars Test 2: returns TRUE if input were created correctly", {
  expect_true(is_order_vars(exprs(AVAL, desc(ADT))))
})

# is_valid_dtc ----
## Test 3: returns TRUE if input are valid dtc ----
test_that("is_valid_dtc Test 3: returns TRUE if input are valid dtc", {
  expect_true(is_valid_dtc("2020"))
  expect_true(is_valid_dtc("2022-09"))
  expect_true(is_valid_dtc("2021-04-06"))
  expect_true(is_valid_dtc("2003-12-15T13:15"))
  expect_true(is_valid_dtc("2021-03-09T01:20:30"))
})

# is_valid_dtc ----
## Test 4: returns error if input if input are NOT valid dtc  ----
test_that("is_valid_dtc Test 4: returns error if input if input are NOT valid dtc ", {
  expect_false(is_valid_dtc("2021-03-T01:20:30"))
})
