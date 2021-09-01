# Re-implemented functions from {lifecycle} package for internal use in case
# the package is not available

deprecate_warn <- function(when, what, with, ...) {
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(when, what, with, ...)
  } else {
    msg <- sprintf(
      "`%s` is deprecated as of admiral %s.\nPlease use `%s` instead.",
      what, when, with
    )
    warn(msg)
  }
}

deprecated <- function() {
  rlang::missing_arg()
}

is_present <- function(arg) {
  !rlang::is_missing(rlang::maybe_missing(arg))
}
