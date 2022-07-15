# ---- my_first_fcn, test 1: hello admiral greets without hw ----
test_that("my_first_fcn, test 1: hello admiral greets without hw", {
  expect_message(
    hello_admiral(),
    "^Welcome to the admiral family!\\n"
  )
})

# ---- my_first_fcn, test 2: hello admiral greets with hw ----
test_that("my_first_fcn, test 2: hello admiral greets with hw", {
  expect_message(
    hello_admiral(TRUE),
    "^Welcome to the admiral family!\\n"
  )
})
