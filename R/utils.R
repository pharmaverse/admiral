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
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' vars2chr(vars(USUBJID, AVAL))
vars2chr <- function(quosures) {
  rlang::set_names(
    map_chr(quosures, ~as_string(quo_get_expr(.x))),
    names(quosures)
  )
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
#' @author Thomas Neitmann, Ondrej Slama
#'
#' @keywords dev_utility
#'
#' @examples
#' test_fun <- function(something) {
#'   admiral:::arg_name(substitute(something))
#' }
#'
#' inner_function <- function(x) x
#' test_fun2 <- function(something) {
#'   admiral:::arg_name(substitute(inner_function(something)))
#' }
arg_name <- function(expr) {
  if (length(expr) == 1L && is.symbol(expr)) {
    deparse(expr)
  } else if (length(expr) == 2L &&
             (expr[[1L]] == quote(enquo) || expr[[1L]] == quote(rlang::enquo)) &&
             is.symbol(expr[[2L]])) {
    deparse(expr[[2L]])
  } else if (is.call(expr) && length(expr) >= 2 && is.symbol(expr[[2]])) {
    deparse(expr[[2L]])
  } else if (is.call(expr) && length(expr) >= 2 && is.call(expr[[2]])) {
    arg_name(expr[[2L]])
  } else {
    abort(paste0("Could not extract argument name from `", deparse(expr), "`"))
  }
}

#' Extract All Symbols from a List of Quosures
#'
#' @param x An `R` object
#' @param side One of `"lhs"` (the default) or `"rhs"`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::extract_vars(vars(STUDYID, USUBJID, desc(ADTM)))
extract_vars <- function(x, side = "lhs") {
  if (is.list(x)) {
    do.call(quo_c, map(x, extract_vars, side))
  } else if (is_quosure(x)) {
    env <- quo_get_env(x)
    symbols <- syms(all.vars(quo_get_expr(x)))
    map(symbols, ~quo_set_env(quo(!!.x), env))
  } else  if (is_formula(x)) {
    funs <- list("lhs" = f_lhs, "rhs" = f_rhs)
    assert_character_scalar(side, values = names(funs))
    quo_set_env(
      quo(!!funs[[side]](x)),
      env = attr(x, ".Environment")
    )
  } else {
    abort()
  }
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

#' Negate List of Variables
#'
#' The function adds a minus sign as prefix to each variable.
#'
#' This is useful if a list of variables should be removed from a dataset,
#' e.g., `select(!!!negate_vars(by_vars))` removes all by variables.
#'
#' @param vars List of variables created by `vars()`
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' negate_vars(vars(USUBJID, STUDYID))
negate_vars <- function(vars = NULL) {
  assert_vars(vars, optional = TRUE)
  if (is.null(vars)) {
    NULL
  } else {
    lapply(vars, function(var) expr(-!!quo_get_expr(var)))
  }
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
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' library(admiral.test)
#' data(vs)
#'
#' filter_if(vs, rlang::quo(NULL))
#' filter_if(vs, rlang::quo(VSTESTCD == "Weight"))
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset)
  assert_filter_cond(filter, optional = TRUE)

  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}

#' Get constant variables
#'
#' @param dataset A data frame.
#' @param by_vars By variables
#'   The groups defined by the by variables are considered separately. I.e., if
#'   a variable is constant within each by group, it is returned.
#'
#' @param ignore_vars Variables to ignore
#'   The specified variables are not considered, i.e., they are not returned
#'   even if they are constant (unless they are included in the by variables).
#'
#'   *Permitted Values:* A list of variable names or selector function calls
#'   like `starts_with("EX")`
#'
#' @keywords dev_utility
#'
#' @return Variable vector.
#'
#' @examples
#' library(admiral.test)
#' data(vs)
#'
#' admiral:::get_constant_vars(vs, by_vars = vars(USUBJID, VSTESTCD))
#'
#' admiral:::get_constant_vars(
#'   vs,
#'   by_vars = vars(USUBJID, VSTESTCD),
#'   ignore_vars = vars(DOMAIN, tidyselect::starts_with("VS"))
#' )
get_constant_vars <- function(dataset, by_vars, ignore_vars = NULL) {
  non_by_vars <- setdiff(names(dataset), vars2chr(by_vars))

  if (!is.null(ignore_vars)) {
    non_by_vars <- setdiff(non_by_vars,
                           vars_select(non_by_vars, !!!ignore_vars))
  }

  # get unique values within each group by variables
  unique_count <- dataset %>%
    group_by(!!!by_vars) %>%
    summarise_at(vars(!!non_by_vars), n_distinct) %>%
    ungroup() %>%
    select(!!!syms(non_by_vars))

  # determine variables which are constant within each by group
  constant_vars <- unique_count %>%
    map_lgl(~ all(.x == 1)) %>%
    which() %>%
    names() %>%
    syms()

  vars(!!!by_vars, !!!constant_vars)
}

`%or%` <- function(lhs, rhs) {
  tryCatch(lhs, error = function(e) rhs)
}

is_named <- function(x) {
  !is.null(names(x)) && all(names(x) != "")
}

replace_values_by_names <- function(quosures) {
  vars <- map2(quosures, names(quosures), function(q, n) {
    if (n == "") {
      return(q)
    }
    quo_set_env(
      quo(!!as.symbol(n)),
      quo_get_env(q)
    )
  })
  structure(vars, class = "quosures", names = NULL)
}

get_duplicates <- function(x) {
  unique(x[duplicated(x)])
}

#' Extract Unit From Parameter Description
#'
#' Extract the unit of a parameter from a description like "Param (unit)".
#'
#' @param x A parameter description
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' extract_unit("Height (cm)")
#'
#' extract_unit("Diastolic Blood Pressure (mmHg)")
extract_unit <- function(x) {
  assert_character_vector(x)

  x %>%
    str_extract("\\(.+\\)") %>%
    str_remove_all("\\(|\\)")
}

#' Convert Blank Strings Into NAs
#'
#' Turn SAS blank strings into proper R `NA`s.
#'
#' @param x Any R object
#'
#' @details
#' The default methods simply returns its input unchanged. The `character` method
#' turns every instance of `""` into `NA_character_` while preserving *all* attributes.
#' When given a data frame as input the function keeps all non-character columns
#' as is and applies the just described logic to `character` columns. Once again
#' all attributes such as labels are preserved.
#'
#' @author Thomas Neitmann
#'
#' @keywords user_utility
#'
#' @export
#'
#' @examples
#' convert_blanks_to_na(c("a", "b", "", "d", ""))
#'
#' df <- tibble::tibble(
#'   a = structure(c("a", "b", "", "c"), label = "A"),
#'   b = structure(c(1, NA, 21, 9), label = "B"),
#'   c = structure(c(TRUE, FALSE, TRUE, TRUE), label = "C"),
#'   d = structure(c("", "", "s", "q"), label = "D")
#' )
#' print(df)
#' convert_blanks_to_na(df)
convert_blanks_to_na <- function(x) {
  UseMethod("convert_blanks_to_na")
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.default <- function(x) {
  x
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.character <- function(x) {
  do.call(structure, c(list(if_else(x == "", NA_character_, x)), attributes(x)))
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.list <- function(x) {
  lapply(x, convert_blanks_to_na)
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.data.frame <- function(x) { # nolint
  x[] <- lapply(x, convert_blanks_to_na)
  x
}

valid_time_units <- function() {
  c("years", "months", "days", "hours", "minutes", "seconds")
}
