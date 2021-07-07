test_that("atomic vectors of length 1", {
  expect_identical(what_is_it(NULL), "`NULL`")
  expect_identical(what_is_it(TRUE), "`TRUE`")
  expect_identical(what_is_it(NA), "`NA`")
  expect_identical(what_is_it("Text"), '`"Text"`')
  expect_identical(what_is_it("3"), '`"3"`')
  expect_identical(what_is_it(4L), "`4`")
  expect_identical(what_is_it(2.42), "`2.42`")
})

test_that("vectors", {
  expect_identical(what_is_it(letters), "a character vector")
  expect_identical(what_is_it(1:10), "an integer vector")
  expect_identical(what_is_it(c(1.2, 3)), "a double vector")
  expect_identical(what_is_it(c(TRUE, FALSE)), "a logical vector")
  expect_identical(what_is_it(list(1, letters, TRUE)), "a list")
})

test_that("S3 objects", {
  expect_identical(what_is_it(mtcars), "a data frame")
  expect_identical(what_is_it(factor(letters)), "a factor")
  expect_identical(what_is_it(lm(hp ~ mpg, data = mtcars)), "an object of class 'lm'")
  expect_identical(what_is_it(quo(4 / 1)), "an object of class 'quosure'")
})


test_that("S4 objects", {
  expect_identical(what_is_it(lubridate::days(1)), "a S4 object of class 'Period'")
})

test_that("symbols", {
  expect_identical(what_is_it(quote(USUBJID)), "a symbol")
})
