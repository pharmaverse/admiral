# get_dataset ----
## Test 1: get_dataset works ----
test_that("get_dataset Test 1: get_dataset works", {
  expect_equal(NULL, get_dataset("one_to_many"))
})

## Test 2: get_dataset gives error  works ----
test_that("get_dataset Test 2: get_dataset works", {
  expect_error(
    get_dataset("test"),
    "`name` must be one of 'one_to_many' or 'many_to_one' but is 'test'"
  )
})
