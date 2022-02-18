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
#' @noRd
#'
#' @examples
#' refdate <- lubridate::ymd("2020-01-02")
#' date <- lubridate::ymd("2020-02-03")
#' assertthat::assert_that(admiral:::is_date(refdate), admiral:::is_date(date))
is_date <- function(arg) {
  is.instant(arg)
}
on_failure(is_date) <- function(call, env) {
  evld <- eval(call$arg, envir = env)
  len <- length(evld)
  msg <- if (len == 0) {
    deparse(evld)
  } else if (len == 1) {
    evld
  } else {
    paste0("c(", paste(head(evld, 5), collapse = ", "), `if`(len > 5, ", ..."), ")")
  }
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    msg,
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
#' @noRd
#'
#' @examples
#' unit <- "days"
#' assertthat::assert_that(admiral:::is_timeunit(unit))
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

#' Check Validity of the Date Imputation Input
#'
#' Date_imputation format should be specified as "dd-mm" (e.g. "01-01")
#' or as a keyword: "FIRST", "MID", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid date_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_date_entry("01-02"))
#' assertthat::assert_that(admiral:::is_valid_date_entry("FIRST"))
is_valid_date_entry <- function(arg) {
  pattern <- "^(01|02|03|04|05|06|07|08|09|10|11|12)-([0-9]{2})$"
  grepl(pattern, arg) | str_to_upper(arg) %in% c("FIRST", "MID", "LAST")
}
on_failure(is_valid_date_entry) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid date entry.\n",
    "date_imputation should be specified as 'mm-dd' (e.g. '01-21') or ",
    "'FIRST', 'MID', 'LAST' to get the first/mid/last day/month"
  )
}

#' Check Validity of the Time Imputation Input
#'
#' Time_imputation format should be specified as "hh:mm:ss" (e.g. "00:00:00")
#' or as a keyword: "FIRST", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid time_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_time_entry("23:59:59"))
#' assertthat::assert_that(admiral:::is_valid_time_entry("FIRST"))
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

#' Check Validity of the Minute/Second Portion of the Time Input
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
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_sec_min(59))
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

#' Check Validity of the Hour Portion in the Time Input
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
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_hour(20))
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

#' Check Validity of the Day Portion in the Date Input
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
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_day(20))
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

#' Check Validity of the Month Portion in the Date Input
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
#' @noRd
#'
#' @examples
#' assertthat::assert_that(admiral:::is_valid_month(12))
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

#' Check Whether an Argument Is Not a Quosure of a Missing Argument
#'
#' @param x Test object
#'
#' @return TRUE or error.
#'
#' @author Thomas Neitmann, Ondrej Slama
#'
#' @noRd
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
