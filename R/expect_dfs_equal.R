#' Expectation: Are Two Datasets Equal?
#'
#' Uses [diffdf::diffdf()] to compares 2 datasets for any differences. This function can be
#' thought of as an R-equivalent of SAS proc compare and a useful tool for unit testing as well.
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
#' @keywords test_helper
#' @family test_helper
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tibble)
#'
#' tbl1 <- tribble(
#'   ~USUBJID, ~AGE, ~SEX,
#'   "1001", 18, "M",
#'   "1002", 19, "F",
#'   "1003", 20, "M",
#'   "1004", 18, "F"
#' )
#'
#' tbl2 <- tribble(
#'   ~USUBJID, ~AGE, ~SEX,
#'   "1001", 18, "M",
#'   "1002", 18.9, "F",
#'   "1003", 20, NA
#' )
#'
#' try(expect_dfs_equal(tbl1, tbl2, keys = "USUBJID"))
#'
#' tlb3 <- tribble(
#'   ~USUBJID, ~AGE, ~SEX,
#'   "1004", 18, "F",
#'   "1003", 20, "M",
#'   "1002", 19, "F",
#'   "1001", 18, "M",
#' )
#'
#' # Note the sorting order of the keys is not required
#' expect_dfs_equal(tbl1, tlb3, keys = "USUBJID")
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
