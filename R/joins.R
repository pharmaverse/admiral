#' Join Functions
#'
#' The `*_join()` functions from {dplyr} without a warning on different attributes
#' in datasets.
#'
#' @inheritParams dplyr::inner_join
#'
#' @return `data.frame`
#'
#' @rdname joins
#' @export
anti_join <- function(x, y, by = NULL, copy = FALSE, ...) {
  suppress_warning(
    dplyr::anti_join(x, y, by = by, copy = copy, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}

#' @rdname joins
#' @export
inner_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::inner_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}

#' @rdname joins
#' @export
left_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}
