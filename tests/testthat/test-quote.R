test_that("Test 1: enumerate works", {
  expect_equal(enumerate(letters[1]), "`a`")
  expect_equal(enumerate(letters[1:3]), "`a`, `b` and `c`")
})

test_that("Test 2: squote works", {
  expect_equal(squote(letters[1]), "'a'")
  expect_equal(squote(letters[1:3]), c("'a'", "'b'", "'c'"))
})

test_that("Test 3: squote works", {
  expect_equal(dquote(letters[1]), "\"a\"")
  expect_equal(dquote(letters[1:3]), c("\"a\"", "\"b\"", "\"c\""))
  x <- NULL
  expect_equal(dquote(x), "NULL")
})
