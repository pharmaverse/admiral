#' Does a Dataset Contain All Required Variables?
#'
#' Checks if a dataset contains all required variables
#'
#' @param dataset A `data.frame`
#' @param required_vars A `character` vector of variable names
#'
#' @author Thomas Neitmann
#'
#' @return The function throws an error if any of the required variables are
#' missing in the input dataset
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(dm)
#' assert_has_variables(dm, "STUDYID")
#' \dontrun{
#' assert_has_variables(dm, "AVAL")
#' }
assert_has_variables <- function(dataset, required_vars) {
  is_missing <- !required_vars %in% colnames(dataset)
  if (any(is_missing)) {
    missing_vars <- required_vars[is_missing]
    if (length(missing_vars) == 1L) {
      err_msg <- paste0("Required variable `", missing_vars, "` is missing.")
    } else {
      err_msg <- paste0(
        "Required variables ",
        enumerate(missing_vars),
        " are missing."
      )
    }
    abort(err_msg)
  }
}

#' Is Date/Date-time?
#'
#' Checks if a date or date-time vector was specified
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a date or date-time, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' refdate <- lubridate::ymd("2020-01-02")
#' date <- lubridate::ymd("2020-02-03")
#' assertthat::assert_that(is_date(refdate), is_date(date))
is_date <- function(arg) {
  is.instant(arg)
}
on_failure(is_date) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a lubridate date."
  )
}

#' Is Time Unit?
#'
#' Checks if a string is a time unit, i.e., 'years', 'months', 'days', 'hours',
#' 'minutes', or 'seconds'.
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a time unit, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' unit <- "days"
#' assertthat::assert_that(is_timeunit(unit))
is_timeunit <- function(arg) {
  arg %in% c("years", "months", "days", "hours", "minutes", "seconds")
}
on_failure(is_timeunit) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid time unit.",
    " Valid time units are 'years', 'months', 'days', 'hours', 'minutes', and 'seconds'."
  )
}

#' Check validity of the date imputation input
#'
#' Date_imputation format should be specified as "dd-mm" (e.g. "01-01")
#' or as a keyword: "FISRT", "MID", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid date_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_date_entry("01-02"))
#' assertthat::assert_that(is_valid_date_entry("FIRST"))
is_valid_date_entry <- function(arg) {
  pattern <- "^([0-9]{2})-([0-9]{2})$"
  grepl(pattern, arg) | str_to_upper(arg) %in% c("FIRST", "MID", "LAST")
}
on_failure(is_valid_date_entry) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid date entry.\n",
    "date_imputation should be specified as 'dd-mm' (e.g. '01-01') or ",
    "'FIRST', 'MID', 'LAST' to get the first/mid/last day/month"
  )
}

#' Check validity of the time imputation input
#'
#' Time_imputation format should be specified as "hh:mm:ss" (e.g. "00:00:00")
#' or as a keyword: "FISRT", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid time_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_time_entry("23:59:59"))
#' assertthat::assert_that(is_valid_time_entry("FIRST"))
is_valid_time_entry <- function(arg) {
  pattern <- "^([0-9]{2}):([0-9]{2}):([0-9]{2})$"
  grepl(pattern, arg) | str_to_upper(arg) %in% c("FIRST", "LAST")
}
on_failure(is_valid_time_entry) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid time entry.\n",
    "time_imputation should be specified as 'hh:mm:ss' (e.g. '00:00:00') or ",
    "'FIRST', 'LAST' to get the first/last time of the day"
  )
}

#' Check validity of the minute/second portion in the time input
#'
#' Minutes and seconds are expected to range from 0 to 59
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid min/sec input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_sec_min(59))
is_valid_sec_min <- function(arg) {
  arg %in% 0:59
}
on_failure(is_valid_sec_min) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid min/sec.\n",
    "Values must be between between 0-59"
  )
}

#' Check validity of the hour portion in the time input
#'
#' Hours are expected to range from 0 to 23
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid hour input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_hour(20))
is_valid_hour <- function(arg) {
  arg %in% 0:23
}
on_failure(is_valid_hour) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    "=",
    eval(call$arg, envir = env),
    " is not a valid hour.\n",
    "Values must be between 0-23"
  )
}

#' Check validity of the day portion in the date input
#'
#' Days are expected to range from 1 to 31
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a day input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_day(20))
is_valid_day <- function(arg) {
  arg %in% 1:31
}
on_failure(is_valid_day) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid day.\n",
    "Values must be between 1-31"
  )
}

#' Check validity of the month portion in the date input
#'
#' Days are expected to range from 1 to 12
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a month input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_month(12))
is_valid_month <- function(arg) {
  arg %in% 1:12
}
on_failure(is_valid_month) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid month.\n",
    "Values for month must be between 1-12. ",
    "Please check the date_imputation input: it should be sepcified as 'dd-mm'"
  )
}

is_named_exprs <- function(arg) {
  is.list(arg) &&
    all(map_lgl(arg, is.language)) &&
    all(names(arg) != "")
}
on_failure(is_named_exprs) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not a named list of expressions created using `exprs()`"
  )
}

#' Is Variable-value List?
#'
#' Checks if the argument is a list of quosures where the expressions are
#' variable-value pairs. The value can be a symbol, a string, or NA. More general
#' expression are not allowed.
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a variable-value list, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_varval_list(vars(DTHDOM = "AE", DTHSEQ = AESEQ)))
is_varval_list <- function(arg) {
  if (inherits(arg, "quosures") && all(names(arg) != "")) {
    expr_list <- map(arg, quo_get_expr)
    all(map_lgl(expr_list, function(arg) is.symbol(arg) || is.character(arg) || is.na(arg)))
  }
  else {
    FALSE
  }
}
on_failure(is_varval_list) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " is not a variable-value pairs list.\n",
    "A named list of quosures is expected where the expression is ",
    "a symbol, a character, or `NA`.\n",
    "The following was supplied:\n",
    paste(capture.output(print(eval(call$arg, envir = env))), collapse = "\n")
  )
}


is_vars <- function(arg) {
  inherits(arg, "quosures") && all(map_lgl(arg, quo_is_symbol))
}
on_failure(is_vars) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not a list of variables created using `vars()`"
  )
}

is_order_vars <- function(arg) {
  quo_is_desc_call <- function(quo) {
    expr <- quo_get_expr(quo)
    is_call(expr) &&
      length(expr) == 2L &&
      deparse(expr[[1L]]) == "desc" &&
      is_symbol(expr[[2L]])
  }

  inherits(arg, "quosures") &&
    all(map_lgl(arg, ~quo_is_symbol(.x) || quo_is_desc_call(.x)))
}
on_failure(is_order_vars) <- function(call, env) {
  paste0(
    backquote(deparse(call$arg)),
    " is not a valid input for `order_vars`.",
    " Valid inputs are created using `vars()` and may only contain symbols or calls involving `desc()`.\n\n", # nolint
    "  # Bad:\n",
    "  vars(ADT = impute_dtc(LBDTC), is.na(AVAL))\n\n",
    "  # Good:\n",
    "  vars(AVAL, desc(ADT))"
  )
}

is_unnamed_exprs <- function(arg) {
  is.list(arg) &&
    all(map_lgl(arg, is.language)) &&
    all(names(arg) == "")
}
on_failure(is_unnamed_exprs) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not a unnamed list of expressions created using `exprs()`"
  )
}

is_expr <- function(arg) {
  # Note: is.language allows both symbol and language
  !is.list(arg) & is.language(arg)
}
on_failure(is_expr) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not an expression created using `expr()`"
  )
}

#' Checks the length of derived records and new values of are equal
#'
#' @param x an R object
#' @param y an R object to compare the length with `x`.
#' @param x_arg Argument name of x.
#' @param y_arg Argument name of y.
#'
#' @return Logical value.
#'
#' @examples
#' \dontrun{
#' x <- list("x", "y")
#' y <- list("y", "z")
#' assertthat::assert_that(are_records_same(x, y, "x", "y"))
#' }
are_records_same <- function(x, y, x_arg, y_arg) {
  stopifnot(is.vector(x), is.vector(y))
  length(x) == length(y)
}

on_failure(are_records_same) <- function(call, env) {
  str_glue("`{call$x_arg}` must have consistent length to the new derived records
           of `{call$y_arg}` within `by_vars`.")
}

#' Check whether an argument is not a quosure of a missing argument
#'
#' @param x Test object
#'
#' @return TRUE or error.
#'
#' @author Thomas Neitmann, Ondrej Slama
#'
#' @export
#'
#' @examples
#' test_fun <- function(x) {x <- rlang::enquo(x); assertthat::assert_that(quo_not_missing(x))}
#' test_fun(my_variable) # no missing argument -> returns TRUE
#' \dontrun{
#' test_fun() # missing argument -> throws error
#' }
quo_not_missing <- function(x) {
  !rlang::quo_is_missing(x)
}
on_failure(quo_not_missing) <- function(call, env) {
  paste0("Argument `", deparse(call$x), "` is missing, with no default")
}
