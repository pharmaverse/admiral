#' @name test_my_first_fcn
#' @title first without using argument
library(testthat)
my_first_fcn <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("")
  }
}

test_that("hello admiral without hw", {
  expect_message(
    my_first_fcn(),
    ""
  )
})
#' @title second with using the argument equal TRUE

test_that("hello admiral with hw", {
  expect_message(
    my_first_fcn(hw= TRUE),
    "Welcome to the admiral family!\\n"
  )
})

