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

#' Are There Multiple Baseline Records?
#'
#' Checks if a dataset contains multiple baseline records
#'
#' @param dataset A `data.frame`
#' @param by A `character` vector of variable names which uniquely identify a
#' set of records that should only contain a single baseline record
#'
#' @author Thomas Neitmann
#'
#' @return The function throws an error if a subject has multiple baseline
#' records
#'
#' @export
assert_has_only_one_baseline_record <- function(dataset, by) { # nolint
  is_duplicate <- duplicated(select(dataset, !!!syms(by)))
  if (any(is_duplicate)) {
    duplicates <- dataset %>%
      select(!!!syms(by)) %>%
      filter(is_duplicate)
    tbl <- capture.output(print(duplicates))
    err_msg <- paste0(
      "Dataset contains multiple baseline records.\n",
      paste(tbl[-c(1, 3)], collapse = "\n")
    )
    abort(err_msg)
  }
}

#' Are records unique?
#'
#' Checks if the reords of a dateset are unique with respect to the specified
#' list of by variables and order.
#'
#' @param dataset The input dataset to check
#'
#' @param by_vars List of by variables
#'
#' @param order Order of observation
#'   If the parameter is specified, it is checked if the observations are unique
#'   with respect to the by variables and the order. If the check fails, the
#'   order values are written as variables in the output.
#'
#' @param message Error message
#'   The message to be displayed if the check fails.
#'
#' @param message_type Message type
#'   If `'error'` is specified, an error is issued if the check fails. Otherwise
#'   an warning is issued.
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a date or date-time, `FALSE` otherwise
#'
#' @export
#'
#' @examples
#' data(ex)
#' assert_has_unique_records(ex,
#'                           by_vars = exprs(USUBJID) ,
#'                           order = exprs(desc(EXENDTC)))

assert_has_unique_records <- function(dataset,
                                      by_vars = NULL,
                                      order = NULL,
                                      message = NULL,
                                      message_type = "error") {
  # variables used for check
  all_vars <- list()

  # variables formatted for the message
  all_vars_msg <- list()

  # dataset to check (remove grouping)
  data_ext <- ungroup(dataset)

  if (!is.null(by_vars)) {
    all_vars <- by_vars
    all_vars_msg <- by_vars
  }
  if (!is.null(order)) {
    # add order variables to the input dataset
    order_vars <- order
    names(order_vars) <- paste0("ordvar", seq_len(length(order_vars)))
    data_ext <- data_ext %>%
      mutate(!!!order_vars)

    # add order variables to the variables for check
    all_vars <- append(all_vars, syms(names(order_vars)))

    # create list of variables for the message, order variables are displayed
    # as ordvar<n> = <expression for order>, e.g., ordvar1 = desc(VISITNUM)
    all_vars_msg <- append(all_vars_msg, paste(names(order_vars), "=", order_vars))
  }

  # select variables for check
  data_by <- data_ext %>% select(!!!all_vars)

  # check for duplicates
  is_duplicate <- duplicated(data_by) | duplicated(data_by, fromLast = TRUE)
  if (any(is_duplicate)) {
    # filter out duplicate observations of the input dataset
    duplicates <- data_ext %>%
      filter(is_duplicate)

    # create message
    tbl <- capture.output(print(duplicates))
    if (is.null(message)) {
      message <- paste0("Dataset contains multiple records with respect to ",
                       paste(all_vars_msg, collapse = ", "),
                       ".")
    }
    err_msg <- paste0(
      message,
      "\n",
      paste(tbl[-c(1, 3)], collapse = "\n")
    )

    # issue message
    if (message_type == "error") {
      abort(err_msg)
    } else {
      warn(err_msg)
    }
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
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_date_entry("01-02"))
#' assertthat::assert_that(is_valid_date_entry("FIRST"))
is_valid_date_entry <- function(arg) {
  pattern <- "^([0-9]{2})-([0-9]{2})$"
  grepl(pattern, arg) | arg %in% c("FIRST", "MID", "LAST")
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
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_time_entry("23:59:59"))
#' assertthat::assert_that(is_valid_time_entry("FIRST"))
is_valid_time_entry <- function(arg) {
  pattern <- "^([0-9]{2}):([0-9]{2}):([0-9]{2})$"
  grepl(pattern, arg) | arg %in% c("FIRST", "LAST")
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
