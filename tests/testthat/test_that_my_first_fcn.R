test_that("hello admiral greets without hm", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family\\n"
  )
})

test_that("hello admiral greets without hm", {
  expect_message(
    hello_admiral(hw = TRUE),
    "^Welcome to the admiral family!\\n"
  )
})

test_that("hello admiral greets without hm", {
  expect_message(
    hello_admiral(hw = FALSE),
    "^Welcome to the admiral family\\n"
  )
})
