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
#' \dontrun{assert_has_variables(dm, "AVAL")}
assert_has_variables <- function(dataset, required_vars) {
  is_missing <- !required_vars %in% colnames(dataset)
  if (any(is_missing)) {
    missing_vars <- required_vars[is_missing]
    if (length(missing_vars) == 1L) {
      err_msg <- paste0("Required variable `", missing_vars, "` is missing.")
    } else {
      err_msg <- paste0("Required variables ",
                        enumerate(missing_vars),
                        " are missing.")
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
assert_has_only_one_baseline_record <- function(dataset, by) {
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
#' refdate <- lubridate::ymd('2020-01-02')
#' date <- lubridate::ymd('2020-02-03')
#' assertthat::assert_that(is_date(refdate), is_date(date))
is_date <- function(arg) {
  is.instant(arg)
}

on_failure(is_date) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    "=",
    eval(call$arg, envir = env),
    " is not a lubridate date."
  )
}

#' Is Time Unit?
#'
#' Checks if a string is a time unit, i.e., 'years', 'months', 'days', 'hours', 'minutes', or 'seconds'.
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
#' unit <- 'days'
#' assertthat::assert_that(is_timeunit(unit))
is_timeunit <- function(arg) {
  arg %in% c("years", "months", "days", "hours", "minutes", "seconds")
}

on_failure(is_timeunit) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    "=",
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
#'
#' assertthat::assert_that(is_valid_date_entry("01-02"))
#' assertthat::assert_that(is_valid_date_entry(FIRST"))

is_valid_date_entry<-function(arg){
  pattern<-"^([0-9]{2})-([0-9]{2})$"
  grepl(pattern,arg) | arg %in% c("FIRST", "MID", "LAST")
}
assertthat::on_failure(is_valid_date_entry) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid date entry.\n",
         "Date_imputation should be specified as 'dd-mm' (e.g. '01-01') or 'FIRST', 'MID', 'LAST' to get the first/mid/last day/month")
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
#'
#' assertthat::assert_that(is_valid_time_entry("23:59:59"))
#' assertthat::assert_that(is_valid_time_entry(FIRST"))
#'
is_valid_time_entry<-function(arg){
  pattern<-"^([0-9]{2}):([0-9]{2}):([0-9]{2})$"
  grepl(pattern,arg) | arg %in% c("FIRST", "LAST")
}
assertthat::on_failure(is_valid_time_entry) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid time entry.\n",
         "time_imputation should be specified as 'hh:mm:ss' (e.g. '00:00:00') or 'FIRST','LAST' to get the first/last time of teh day")
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
#'
#' assertthat::assert_that(is_valid_sec_min(59))
#'
is_valid_sec_min<-function(arg){
  arg %in% 0:59
}
assertthat::on_failure(is_valid_sec_min) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid min/sec.\n",
         "Values must be between between 0-59")
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
#'
#' assertthat::assert_that(is_valid_hour(20))
#'

is_valid_hour<-function(arg){
  arg %in% 0:23
}
assertthat::on_failure(is_valid_hour) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid hour.\n",
         "Values must be between 0-23")
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
#'
#' assertthat::assert_that(is_valid_day(20))
#'

is_valid_day<-function(arg){
  arg %in% 1:31
}
assertthat::on_failure(is_valid_day) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid day.\n",
         "Values must be between 1-31")
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
#'
#' assertthat::assert_that(is_valid_month(20))
#'

is_valid_month<-function(arg){
  arg %in% 1:12
}
assertthat::on_failure(is_valid_month) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid month.\n",
         "Values for month must be between 1-12. Please check the date_imputation input: it should be sepcified as 'dd-mm'")
}
