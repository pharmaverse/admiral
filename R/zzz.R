
# Re-implemented functions from "lifecycle" package

deprecate_warn <- function(...) {
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(...)
  } else {
    invisible()
  }
}

deprecated <- function() rlang::missing_arg()

is_present <- function(arg) !rlang::is_missing(rlang::maybe_missing(arg))
