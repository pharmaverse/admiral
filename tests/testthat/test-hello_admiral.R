#  hello_admiral -------
## Test 1: hello_admiral works ----
test_that("hello_admiral() Test #1  works", {
  expect_no_error(hello_admiral(TRUE))
  expect_no_error(hello_admiral(FALSE))
  expect_no_error(hello_admiral(0))
  expect_no_error(hello_admiral())
})

## Test2: hello_admiral  throws Error if input not acceptable ----
test_that("hello_admiral() Test #2 Error if input not acceptable type", {
  expect_error(hello_admiral(""))
  expect_error(hello_admiral("junk"))
  expect_error(hello_admiral(junk))
})

## Test3: hello_admiral ..... ----
test_that("hello_admiral() Test #3: Return value is `interpreted string literal`", {
  s <- "Welcome to Admiral family"
  expect_identical(hello_admiral(), cat(s, "\n"))
})
