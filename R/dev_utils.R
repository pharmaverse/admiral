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
#' @rdname dev_util_enumerate
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
#' @rdname dev_util_backquote
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
#' @rdname dev_util_squote
#'
#' @examples
#' admiral:::squote("foo")
squote <- function(x) {
  paste0("'", x, "'")
}

#' Wrap a String in Double Quotes
#'
#' Wrap a string in double quotes, e.g., for displaying character values in
#' messages.
#'
#' @param x A character vector
#'
#' @return If the input is `NULL`, the text `"NULL"` is returned. Otherwise, the
#'   input in double quotes is returned.
#'
#' @author Stefan Bundfuss
#'
#' @keywords dev_utility
#'
#' @examples
#' admiral:::dquote("foo")
#' admiral:::dquote(NULL)
dquote <- function(x) {
  if (is.null(x)) {
    "NULL"
  }
  else {
    paste0("\"", x, "\"")
  }
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
#' @rdname dev_util_notin
#'
#' @examples
#' `%notin%` <- admiral:::`%notin%`
#' "a" %notin% c("b", "v", "k")
`%notin%` <- function(x, table) { # nolint
  !(x %in% table)
}

#' Helper Function to Convert Date (or Date-time) Objects to Characters of dtc Format
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
#' @rdname dev_util_convert_dtm_to_dtc
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
#' @rdname dev_util_arg_name
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
arg_name <- function(expr) { # nolint
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
#' @return A list of `quosures`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_extract_vars
#'
#' @examples
#' admiral:::extract_vars(vars(STUDYID, USUBJID, desc(ADTM)))
extract_vars <- function(x, side = "lhs") {
  if (is.null(x)) {
    NULL
  } else if (is.list(x)) {
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
#' @return An object of class `quosures`
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_quo_c
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
#' @rdname dev_util_what_is_it
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
    friendly_type_of(x)
  } else if (is.data.frame(x)) {
    "a data frame"
  } else {
    sprintf("an object of class '%s'", class(x)[1L])
  }
}

#' Get Constant Variables
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
#' @rdname dev_util_get_constant_vars
#'
#' @return Variable vector.
#'
#' @examples
#' library(admiraltest)
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

#' Replace Quosure Value with Name
#'
#' @param quosures A list of quosures
#'
#' @author Thomas Neitmann
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_replace_values_by_names
#'
#' @return A list of quosures
#'
#' @examples
#' admiral:::replace_values_by_names(vars(USUBJID, TEST = VSTESTCD))
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

#' Checks if the argument equals the auto keyword
#'
#' @param arg argument to check
#'
#' @return `TRUE` if the argument equals the auto keyword, i.e., it is a quosure
#'   of a symbol named auto.
#'
#' @author Stefan Bundfuss
#'
#' @keywords check
#'
#' @examples
#'
#' example_fun <- function(arg) {
#'   arg <- rlang::enquo(arg)
#'   if (admiral:::is_auto(arg)) {
#'     "auto keyword was specified"
#'   }
#'   else {
#'     arg
#'   }
#' }
#'
#' example_fun("Hello World!")
#'
#' example_fun(auto)
is_auto <- function(arg) {
  is_quosure(arg) && quo_is_symbol(arg) && quo_get_expr(arg) == expr(auto)
}

#' Get Source Variables from a List of Quosures
#'
#' @param quosures A list of quosures
#'
#' @author Stefan Bundfuss
#'
#' @keywords dev_utility
#'
#' @rdname dev_util_get_source_vars
#'
#' @return A list of quosures
#'
#' @examples
#' admiral:::get_source_vars(vars(USUBJID, AVISIT = VISIT, SRCDOM = "EX"))
get_source_vars <- function(quosures) {
  quo_c(quosures)[lapply(quo_c(quosures), quo_is_symbol) == TRUE]
}

valid_time_units <- function() {
  c("years", "months", "days", "hours", "minutes", "seconds")
}

contains_vars <- function(arg) {
  inherits(arg, "quosures") && all(map_lgl(arg, quo_is_symbol) | names(arg) != "")
}

#' Get a New Temporary Variable Name for a Dataset
#'
#' @param dataset The input dataset
#' @param prefix The prefix of the new temporary variable name to create
#'
#' @details
#' The function returns a new unique temporary variable name to be used inside
#' `dataset`. The temporary variable names have the structure `prefix_n` where
#' `n` is an integer, e.g. `tmp_var_1`. If there is alredy a variable inside
#' `datset` with a given `prefix` then the suffix is increased by 1, e.g. if
#' `tmp_var_1` already exists then `get_new_tmp_var()` will return `tmp_var_2`.
#'
#' @examples
#' data(adsl)
#' tmp_var <- get_new_tmp_var(adsl)
#' mutate(adsl, !!tmp_var := NA)
get_new_tmp_var <- function(dataset, prefix = "tmp_var") {
  assert_data_frame(dataset)
  assert_character_scalar(prefix)
  if (!str_detect(prefix,"^tmp_")) {
    abort("`prefix` must start with 'tmp_'")
  }

  caller_env <- parent.frame()

  regexp <- str_c("^", prefix, "_[0-9]{1,}$")
  colnames_with_prefix <- str_subset(colnames(dataset), regexp)
  if (!is.null(caller_env$.tmp_vars)) {
    suffices <- str_extract(caller_env$.tmp_vars, "[0-9]{1,}$")
    counter <- max(as.integer(suffices)) + 1L
  } else if (length(colnames_with_prefix) > 0L) {
    suffices <- str_extract(colnames_with_prefix, "[0-9]{1,}$")
    counter <- max(as.integer(suffices)) + 1L
  } else {
    counter <- 1L
  }

  new_tmp_var <- str_c(prefix, counter, sep = "_")
  caller_env$.tmp_vars <- c(caller_env$.tmp_vars, new_tmp_var)

  sym(new_tmp_var)
}

#' Remove All Temporary Variables Created Within the Current Function Environment
#'
#' @param dataset The input dataset
#'
#' @examples
#' data(dm)
#' tmp_var <- get_new_tmp_var(dm)
#' dm <- mutate(dm, !!tmp_var := NA)
#'
#' ## This function creates two new temporary variables which are removed when calling
#' ## `remove_tmp_vars()`. Note that any temporary varirable created outside this
#' ## function is **not** removed
#' do_something <- function(dataset) {
#'   tmp_var_1 <- get_new_tmp_var(dm)
#'   tmp_var_2 <- get_new_tmp_var(dm)
#'   dm %>%
#'     mutate(!!tmp_var_1 := NA, !!tmp_var_2 := NA) %>%
#'     print() %>%
#'     remove_tmp_vars()
#' }
#'
#' do_something(dm)
remove_tmp_vars <- function(dataset) {
  # In order to find the "correct" calling environment we have to make sure to
  # exclude all calls relating to the use of `%>%` from the call stack
  calls <- lapply(sys.calls(), function(x) paste(deparse(x), collapse = " "))
  contains_pipe <- map_lgl(calls, str_detect, "%>%")

  if (any(contains_pipe)) {
    last_pipe <- max(which(contains_pipe))
    n <- length(calls) - last_pipe - 1L
  } else {
    n <- 1L
  }

  tmp_vars <- parent.frame(n)$.tmp_vars
  if (is.null(tmp_vars)) {
    dataset
  } else {
    dataset[, colnames(dataset) %notin% tmp_vars]
  }
}
