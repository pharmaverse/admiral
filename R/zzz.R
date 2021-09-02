
# Re-implemented functions from "lifecycle" package

deprecate_warn <- function(...) {
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(...)
  } else {
    invisible()
  }
}

deprecated <- function() rlang::missing_arg()
