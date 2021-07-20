#' Enumerate Multiple Strings
#'
#' @param x A `character` vector
#' @param quote_fun Quoting function, defaults to `backquote`.
#' @param conjunction Character to be used in the message, defaults to "and".
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::enumerate(c("STUDYID", "USUBJID", "PARAMCD"))
#' admiral:::enumerate(letters[1:6], quote_fun = admiral:::squote)
#' admiral:::enumerate(
#'   c("date", "time", "both"),
#'   quote_fun = admiral:::squote,
#'   conjunction = "or"
#' )
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
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::backquote("USUBJID")
backquote <- function(x) {
  paste0("`", x, "`")
}

#' Wrap a String in Single Quotes
#'
#' @param x A `character` vector
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::squote("foo")
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
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' `%notin%` <- admiral:::`%notin%`
#' "a" %notin% c("b", "v", "k")
`%notin%` <- function(x, table) { # nolint
  !(x %in% table)
}

#' Turn a List of Quosures into a Character Vector
#'
#' @param quosures A `list` of `quosures` created using [`vars()`]
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::vars2chr(vars(USUBJID, AVAL))
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
#' @author Ondrej Slama
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::convert_dtm_to_dtc(as.POSIXct(Sys.time()))
#' admiral:::convert_dtm_to_dtc(as.Date(Sys.time()))
convert_dtm_to_dtc <- function(dtm) {
  stopifnot(lubridate::is.instant(dtm))
  format(dtm, "%Y-%m-%dT%H:%M:%S")
}

#' Extract Argument Name from an Expression
#'
#' @param expr An expression created inside a function using `substitute()`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' test_fun <- function(something) {
#'   admiral:::arg_name(substitute(something))
#' }
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

#' Extract All Symbols from a List of Quosures
#'
#' @param x An `R` object
#' @param side One of `"lhs"` (the default) or `"rhs"`
#' @param ... Not used
#'
#' @author Thomas Neitmann
#' @export
#' @keywords dev_utility
#'
#' @examples
#' extract_vars(vars(STUDYID, USUBJID, desc(ADTM)))
extract_vars <- function(x, ...) {
  UseMethod("extract_vars")
}

#' @describeIn extract_vars
#' @export
extract_vars.list <- function(x, ...) {
  do.call(quo_c, map(x, extract_vars, ...))
}

#' @describeIn extract_vars
#' @export
extract_vars.quosure <- function(x, ...) {
  env <- quo_get_env(x)
  symbols <- syms(all.vars(quo_get_expr(x)))
  map(symbols, ~quo_set_env(quo(!!.x), env))
}

#' @describeIn extract_vars
#' @export
extract_vars.formula <- function(x, side = "lhs", ...) {
  funs <- list("lhs" = f_lhs, "rhs" = f_rhs)
  assert_character_scalar(side, values = names(funs))
  quo_set_env(
    quo(!!funs[[side]](x)),
    env = attr(x, ".Environment")
  )
}

#' Concatenate One or More Quosure(s)
#'
#' @param ... One or more objects of class `quosure` or `quosures`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::quo_c(rlang::quo(USUBJID))
#' admiral:::quo_c(rlang::quo(STUDYID), rlang::quo(USUBJID))
#' admiral:::quo_c(vars(USUBJID, ADTM))
#' admiral:::quo_c(rlang::quo(BASETYPE), vars(USUBJID, PARAM), rlang::quo(ADTM))
quo_c <- function(...) {
  inputs <- unlist(list(...), recursive = TRUE)
  stopifnot(all(map_lgl(inputs, is_quosure)))
  is_null <- map_lgl(inputs, quo_is_null)
  rlang::as_quosures(inputs[!is_null])
}

#' What Kind of Object is This?
#'
#' Returns a string describing what kind of object the input is.
#'
#' @param x Any R object
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::what_is_it(mtcars)
#' admiral:::what_is_it(NA)
#' admiral:::what_is_it(TRUE)
#' admiral:::what_is_it(lm(hp ~ mpg, data = mtcars))
#' admiral:::what_is_it(letters)
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

#' Optional Filter
#'
#' Filters the input dataset if the provided expression is not `NULL`
#'
#' @param dataset Input dataset
#' @param filter A filter condition. Must be a quosure.
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' data(vs)
#' admiral:::filter_if(vs, rlang::quo(NULL))
#' admiral:::filter_if(vs, rlang::quo(VSTESTCD == "Weight"))
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset)
  assert_filter_cond(filter, optional = TRUE)

  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}
