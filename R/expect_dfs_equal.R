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

#' Expectation: Are Table Lists are Equal?
#'
#' Uses [diffdf::diffdf()] to compares 2 datasets for any differences
#'
#' @param base Input dataset
#' @param compare Comparison dataset

#' @return
#' An error if `base` and `compare` do not match or `NULL` invisibly if they do
#'
#' @author Ben Straub
#' @keywords test_helper
#' @family test_helper
#'
#' @export
expect_equal_tbl <- function(object, expected, ..., info = NULL) {
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
  exp <- testthat::quasi_label(rlang::enquo(expected), arg = "expected")

  # all.equal.list is slightly problematic: it returns TRUE for match, and
  # returns a character vector when differences are observed. We extract
  # both a match-indicator and a failure message

  diffs <- all.equal.list(object, expected, ...)
  has_diff <- if (is.logical(diffs)) diffs else FALSE
  diff_msg <- paste(diffs, collapse = "\n")

  testthat::expect(
    has_diff,
    failure_message = sprintf(
      "%s not equal to %s.\n%s", act$lab, exp$lab, diff_msg
    ),
    info = info
  )

  invisible(act$val)
}
