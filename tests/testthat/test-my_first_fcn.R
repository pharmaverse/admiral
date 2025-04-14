test_that("hello_admiral() Test 1: snapshot test with hw = TRUE", {
  expect_snapshot({
    hello_admiral(hw = TRUE)
  })
})

test_that("hello_admiral() Test 2: snapshot test with hw = FALSE", {
  expect_snapshot({
    hello_admiral(hw = FALSE)
  })
})

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
