#' @name test_my_first_fcn
my_first_fcn <- function(hw = TRUE) {
  if (hw) {
    message("Welcome to the admiral family!")
  } else {
    message("Welcome to the admiral family!")
  }
}
#' @title first without using argument
test_that("hello admiral without hw", {
  expect_message(
    my_first_fcn(),
    "^Welcome to the admiral family!\\n"
  )
})
#' @title second with using the argument equal TRUE

test_that("hello admiral with hw", {
  expect_message(
    my_first_fcn(hw = TRUE),
    "^Welcome to the admiral family!\\n"
  )
})
