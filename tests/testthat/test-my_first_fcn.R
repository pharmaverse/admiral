# hello_admiral ----
## Test 1: hello admiral greets without smile ----
test_that("hello_admiral Test 1: hello admiral greets without smile", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family :)\\n"
  )
})

## Test 2: hello admiral greets with smile ----
test_that("hello_admiral Test 2: hello admiral greets with smile", {
  expect_message(
    hello_admiral(smile = TRUE),
    "^Welcome to the admiral family :)\\n"
  )
})

## Test 3: hello admiral greets with smile ----
test_that("hello_admiral Test 3: hello admiral greets with smile", {
  expect_message(
    hello_admiral(smile = FALSE),
    "^Welcome to the admiral family!\\n"
  )
})
