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
