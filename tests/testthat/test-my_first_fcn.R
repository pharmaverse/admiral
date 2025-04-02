test_that("hello admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "Welcome to the admiral family with hw"
  )
})

test_that("hello admiral greets with hw", {
  expect_message(
    hello_admiral(hw = TRUE),
    "Welcome to the admiral family with hw"
  )
})

test_that("hello admiral greets with hw", {
  expect_message(
    hello_admiral(hw = FALSE),
    "Welcome to the admiral family without hw"
  )
})
