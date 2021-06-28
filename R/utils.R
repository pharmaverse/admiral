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
