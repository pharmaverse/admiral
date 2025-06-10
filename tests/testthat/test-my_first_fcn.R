# hello_admiral ----
## Test 1: hello admiral greets without hw ----
test_that("hello_admiral Test 1: hello admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family!\\n"
  )
})

## Test 2: hello admiral greets with hw ----
test_that("hello_admiral Test 2: hello admiral greets with hw", {
  expect_message(
    hello_admiral(hw = TRUE),
    "^Welcome to the admiral family!\\n"
  )
})

## Test 3: hello admiral greets with hw ----
test_that("hello_admiral Test 3: hello admiral greets with hw", {
  expect_message(
    hello_admiral(hw = FALSE),
    "^Welcome to the admiral family!\\n"
  )
})
