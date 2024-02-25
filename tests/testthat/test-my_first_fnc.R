test_that("hello_admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family!\\n"
  )
})

test_that("hello_admiral greets with hw", {
  expect_message(
    hello_admiral(TRUE),
    "^Welcome to the admiral family!\\n"
  )
})
