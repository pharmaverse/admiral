test_that("hello admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family!\\n"
  )
})

test_that("hello admiral greets with hw", {
  expect_message(
    hello_admiral(TRUE),
    "^Welcome to the admiral family!\\n"
  )
})
