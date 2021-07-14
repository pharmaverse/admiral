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
#' @param quote_fun Quoting function, defaults to `backquote`.
#' @param conjunction Character to be used in the message, defaults to "and".
#'
#' @noRd
#'
#' @examples
#' enumerate(letters[1:6])
enumerate <- function(x, quote_fun = backquote, conjunction = "and") {
  if (length(x) == 1L) {
    quote_fun(x)
  } else {
    paste(
      paste0(quote_fun(x[-length(x)]), collapse = ", "),
      conjunction,
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
#' @noRd
#'
#' @examples
#' "a" %!in% c("b", "v", "k")
`%!in%` <- function(x, table) { # nolint
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

#' Helper function to convert date (or date-time) objects to characters of dtc format
#' (-DTC type of variable)
#'
#' @param dtm date or date-time
#'
#' @return character
#'
#' @noRd
#'
#' @examples
#' admiral:::convert_dtm_to_dtc(as.POSIXct(Sys.time()))
#' admiral:::convert_dtm_to_dtc(as.Date(Sys.time()))
convert_dtm_to_dtc <- function(dtm) {
  stopifnot(lubridate::is.instant(dtm))
  format(dtm, "%Y-%m-%dT%H:%M:%S")
}

arg_name <- function(expr) {
  if (length(expr) == 1L && is.symbol(expr)) {
    deparse(expr)
  } else if (length(expr) == 2L &&
             (expr[[1L]] == quote(enquo) || expr[[1L]] == quote(rlang::enquo)) &&
             is.symbol(expr[[2L]])) {
    deparse(expr[[2L]])
  } else {
    abort(paste0("Could not extract argument name from `", deparse(expr), "`"))
  }
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

quo_c <- function(...) {
  inputs <- unlist(list(...), recursive = TRUE)
  stopifnot(all(map_lgl(inputs, is_quosure)))
  is_null <- map_lgl(inputs, quo_is_null)
  rlang::as_quosures(inputs[!is_null])
}

what_is_it <- function(x) {
  if (is.null(x)) {
    "`NULL`"
  } else if (is.factor(x)) {
    "a factor"
  } else if (is.symbol(x)) {
    "a symbol"
  } else if (isS4(x)) {
    sprintf("a S4 object of class '%s'", class(x)[1L])
  } else if (is.atomic(x) && length(x) == 1L) {
    if (is.character(x)) {
      paste0("`\"", x, "\"`")
    } else {
      paste0("`", x, "`")
    }
  } else if (is.atomic(x) || class(x)[1L] == "list") {
    rlang::friendly_type(typeof(x))
  } else if (is.data.frame(x)) {
    "a data frame"
  } else {
    sprintf("an object of class '%s'", class(x)[1L])
  }
}
