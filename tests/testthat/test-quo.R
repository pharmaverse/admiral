test_that("Test 10 : `quo_not_missing` returns TRUE if no missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_true(test_fun(my_variable))
})

test_that("Test 11 : `quo_not_missing` throws and Error if missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})
