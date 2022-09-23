# enumerate ----
## Test 1: enumerate works ----
test_that("enumerate Test 1: enumerate works", {
  expect_equal(enumerate(letters[1]), "`a`")
  expect_equal(enumerate(letters[1:3]), "`a`, `b` and `c`")
})

# squote ----
## Test 2: squote works ----
test_that("squote Test 2: squote works", {
  expect_equal(squote(letters[1]), "'a'")
  expect_equal(squote(letters[1:3]), c("'a'", "'b'", "'c'"))
})

# dquote ----
## Test 3: dquote works ----
test_that("dquote Test 3: dquote works", {
  expect_equal(dquote(letters[1]), "\"a\"")
  expect_equal(dquote(letters[1:3]), c("\"a\"", "\"b\"", "\"c\""))
  x <- NULL
  expect_equal(dquote(x), "NULL")
})
