#' @export
dplyr::vars

#' @export
dplyr::desc

#' @export
rlang::exprs

#' @export
dplyr::vars

#' Enumerate Multiple Strings
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' enumerate(letters[1:6])
enumerate <- function(x, quote_fun = backquote) {
  if (length(x) == 1L) {
    quote_fun(x)
  } else {
    paste(
      paste0(quote_fun(x[-length(x)]), collapse = ", "),
      "and",
      quote_fun(x[length(x)])
    )
  }
}

#' Wrap a String in Backquotes
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' backquote("foo")
backquote <- function(x) {
  paste0("`", x, "`")
}

#' Wrap a String in Single Quotes
#'
#' @param x A `character` vector
#'
#' @noRd
#'
#' @examples
#' squote("foo")
squote <- function(x) {
  paste0("'", x, "'")
}

#' Negated Value Matching
#'
#' Returns a `logical` vector indicating if there is *no* match of the
#' left operand in the right operand.
#'
#' @param x The values to be matched
#' @param table The values to be matched against
#'
#' @rdname utils
#' @name not_in
#'
#' @export
#'
#' @examples
#' "a" %!in% c("b", "v", "k")
`%!in%` <- function(x, table) {
  !(x %in% table)
}

#' Turn a List of Quosures into a Character Vector
#'
#' @param quosures A `list` of `quosures` created using [`vars()`]
#'
#' @noRd
#'
#' @examples
#' vars2chr(vars(USUBJID, AVAL))
vars2chr <- function(quosures) {
  map_chr(quosures, ~as_string(quo_get_expr(.x)))
}


#' Return y if x has length > 0
#'
#' @param x vector
#' @param y vector
#'
#' @noRd
#'
#' @examples
#' if_non_len0(character(0), "whatever")
#' if_non_len0(character(1), "updated if length > 0")
if_non_len0 <- function(x, y) {
  if (length(x) != 0) {
    y
  } else {
    x
  }
}

#' Helper function to convert date (or date-time) objects to characters of dtc format
#' (-DTC type of variable)
#'
#' @param dtm date or date-time
#'
#' @return character
#'
#' @examples
#' admiral:::convert_dtm_to_dtc(as.POSIXct(Sys.time()))
#' admiral:::convert_dtm_to_dtc(as.Date(Sys.time()))
convert_dtm_to_dtc <- function(dtm) {
  stopifnot(lubridate::is.instant(dtm))
  format(dtm, "%Y-%m-%dT%H:%M:%S")
}

extract_vars <- function(quosures) {
  vars <- lapply(quosures, function(q) {
    rlang::quo_set_env(
      rlang::quo(!!as.symbol(all.vars(q))),
      rlang::quo_get_env(q)
    )
  })
  structure(vars, class = "quosures")
}

left_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}

inner_join <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  suppress_warning(
    dplyr::inner_join(x, y, by = by, copy = copy, suffix = suffix, ...),
    "^Column `.+` has different attributes on LHS and RHS of join$"
  )
}
