
# Re-implemented functions from "lifecycle" package

deprecate_warn <- function(...) {
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(...)
  } else {
    invisible()
  }
}

deprecated <- function() rlang::missing_arg()

is_present <- function(arg, quoted = FALSE) {
  stopifnot(rlang::is_scalar_logical(quoted))
  if (quoted) {
    arg <- rlang::quo_get_expr(arg)
    if (identical(arg, quote(deprecated()))) {
      return(FALSE)
    }
  }
  !rlang::is_missing(rlang::maybe_missing(arg))
}
