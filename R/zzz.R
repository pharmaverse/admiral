deprecate_warn <- function(...) {
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(...)
  } else {
    invisible()
  }
}
