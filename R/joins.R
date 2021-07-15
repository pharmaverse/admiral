inner_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::inner_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}

left_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}
