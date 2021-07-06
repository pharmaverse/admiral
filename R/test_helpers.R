#' Expectation: Are Two Datasets Equal?
#'
#' Uses [diffdf::diffdf()] to compares 2 datasets for any differences
#'
#' @param base Input dataset
#' @param compare Comparison dataset
#' @param keys `character` vector of variables that define a unique row in the
#'        `base` and `compare` datasets
#' @param ... Additional arguments passed onto [diffdf::diffdf()]
#'
#' @author Thomas Neitmann
#' @keywords test_helper
#' @export
#'
#' @examples
#' \dontrun{
#' testthat::test_that("a missing row is detected", {
#'   data(dm)
#'   expect_dfs_equal(dm, dm[-1L, ], keys = "USUBJID")
#' })
#' }
expect_dfs_equal <- function(base, compare, keys, ...) {
  diff <- diffdf::diffdf(base, compare, keys, suppress_warnings = TRUE, ...)
  if (diffdf::diffdf_has_issues(diff)) {
    msg <- capture.output(print(diff))
    testthat::fail(msg)
  } else {
    testthat::succeed()
    invisible()
  }
}
