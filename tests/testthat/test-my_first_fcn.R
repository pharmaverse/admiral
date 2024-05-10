testthat::test_that("hello admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "Welcome to Admiral!"
  )

})


testthat::test_that("hello admiral greets without hw", {
  expect_message(
    hello_admiral(T),
    "Welcome to Admiral!"
  )

})


testthat::test_that("hello admiral greets without hw", {
  expect_message(
    hello_admiral(F),
    "Welcome to Admiral!"
  )

})
