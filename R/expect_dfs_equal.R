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
#' @return
#' An error if `base` and `compare` do not match or `NULL` invisibly if they do
#'
#' @author Thomas Neitmann
#' @keywords test_helper
#' @family test_helper
#'
#' @export
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
