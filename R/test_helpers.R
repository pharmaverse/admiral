#' Expectation: Are Two Datasets Equal
#'
#' Uses [diffdf::diffdf()] to compares 2 datasets for any differences
#'
#' @param base Input dataset
#' @param compare Comparison dataset
#' @param keys `character` vector of variables that define a unique row in the
#'        base and compare datasets
#' @param ... Additional arguments passed onto [diffdf::diffdf()]
#'
#' @examples
#' testthat::test_that("a missing row is detected", {
#'   data(dm)
#'   expect_dfs_equal(dm, dm[-1L, ], keys = "USUBJID")
#' })
#'
expect_dfs_equal <- function(base, compare, keys, ...) {
  diff <- suppressWarnings(diffdf::diffdf(base, compare, keys, ...))
  if (length(diff) == 0L) {
    testthat::succeed()
    invisible()
  } else {
    msg <- capture.output(print(diff))
    testthat::fail(msg)
  }
}
