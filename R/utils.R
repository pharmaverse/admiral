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
  paste(
    paste0(quote_fun(x[-length(x)]), collapse = ", "),
    "and",
    quote_fun(x[length(x)])
  )
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

#' @export
print <- function(x, ...) {
  UseMethod("print")
}

#' @export
print.default <- function(x, ...) {
  base::print(x, ...)
}

#' @export
print.POSIXct <- function (x, tz = "", usetz = FALSE, max = NULL, ...) {
  if (is.null(max)) max <- getOption("max.print", 9999L)

  FORM <- if (missing(tz))
    function(z) format(z, usetz = usetz, format = "%Y-%m-%d %H:%M:%S")
  else
    function(z) format(z, tz = tz, usetz = usetz, format = "%Y-%m-%d %H:%M:%S")

  if (max < length(x)) {
    print(FORM(x[seq_len(max)]), max = max + 1, ...)
    cat(" [ reached 'max' / getOption(\"max.print\") -- omitted",
        length(x) - max, "entries ]\n")
  } else if (length(x)) {
    print(FORM(x), max = max, ...)
  } else {
    cat(class(x)[1L], "of length 0\n")
  }

  invisible(x)
}

