# what_is_it ----
## Test 1: atomic vectors of length 1 ----
test_that("what_is_it Test 1: atomic vectors of length 1", {
  expect_identical(what_is_it(NULL), "`NULL`")
  expect_identical(what_is_it(TRUE), "`TRUE`")
  expect_identical(what_is_it(NA), "`NA`")
  expect_identical(what_is_it("Text"), '`"Text"`')
  expect_identical(what_is_it("3"), '`"3"`')
  expect_identical(what_is_it(4L), "`4`")
  expect_identical(what_is_it(2.42), "`2.42`")
})

## Test 2: vectors ----
test_that("what_is_it Test 2: vectors", {
  expect_identical(what_is_it(letters), "a character vector")
  expect_identical(what_is_it(1:10), "an integer vector")
  expect_identical(what_is_it(c(1.2, 3)), "a double vector")
  expect_identical(what_is_it(c(TRUE, FALSE)), "a logical vector")
  expect_identical(what_is_it(list(1, letters, TRUE)), "a list")
})

## Test 3: S3 objects ----
test_that("what_is_it Test 3: S3 objects", {
  expect_identical(what_is_it(mtcars), "a data frame")
  expect_identical(what_is_it(factor(letters)), "a factor")
  expect_identical(what_is_it(lm(hp ~ mpg, data = mtcars)), "an object of class 'lm'")
  expect_identical(what_is_it(quo(4 / 1)), "an object of class 'quosure'")
})


## Test 4: S4 objects ----
test_that("what_is_it Test 4: S4 objects", {
  expect_identical(what_is_it(lubridate::days(1)), "a S4 object of class 'Period'")
})

## Test 5: symbols ----
test_that("what_is_it Test 5: symbols", {
  expect_identical(what_is_it(quote(USUBJID)), "a symbol")
})
