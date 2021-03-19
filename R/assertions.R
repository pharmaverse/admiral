assert_has_variables <- function(dataset, required_vars) {
  is_missing <- !required_vars %in% colnames(dataset)
  if (any(is_missing)) {
    missing_vars <- required_vars[is_missing]
    if (length(missing_vars) == 1L) {
      err_msg <- paste0("Required variable `", missing_vars, "` is missing.")
    } else {
      err_msg <- paste0("Required variables ", enumerate(missing_vars), " are missing.")
    }
    abort(err_msg)
  }
}

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
#' @return ``TRUE`` if the argument is a date or date-time, ``FALSE`` otherwise
#'
#' @export
#'
#' @examples
#' assert_that(is_date(refdate), is_date(date))

is_date <- function(arg){
  is.instant(arg)
}

assertthat::on_failure(is_date) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a lubridate date.")
}

#' Is Time Unit?
#'
#' Checks if a string is a time unit, i.e., 'years', 'months', 'days', 'hours', 'minutes', or 'seconds'.
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return ``TRUE`` if the argument is a time unit, ``FALSE`` otherwise
#'
#' @export
#'
#' @examples
#' assert_that(is_timeunit(unit))

is_timeunit <- function(arg){
  arg %in% c('years', 'months', 'days', 'hours', 'minutes', 'seconds')
}

assertthat::on_failure(is_timeunit) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid time unit.",
         " Valid time units are 'years', 'months', 'days', 'hours', 'minutes', and 'seconds'.")
}


is_valid_sec_min<-function(arg){
  arg %in% seq(from=0,to=59)
}
is_valid_hour<-function(arg){
  arg %in% seq(from=0,to=23)
}
is_valid_day<-function(arg){
  arg %in% seq(1:31) | arg=="LAST"
}
is_valid_month<-function(arg){
  arg %in% seq(1:12)
}
is_valid_year<-function(arg){
  nchar(arg)==4
}
assertthat::on_failure(is_valid_sec_min) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid min/sec.",
         "Values must be between between 0-59")
}
assertthat::on_failure(is_valid_hour) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid day.",
         "Values must be between 0-23")
}
assertthat::on_failure(is_valid_day) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid day.",
         "Values must be between 1-31 or 'LAST'")
}
assertthat::on_failure(is_valid_month) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid month",
         "Values must be between 1-12")
}
assertthat::on_failure(is_valid_year) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid year",
         "Expecting a 4-digit number")
}

assert_has_a_date_variable <- function(dataset, dates) {
  is_missing <- !dates %in% colnames(dataset)
  if (all(is_missing)) {
    err_msg <- paste("Required a date variable: `", dates[1], dates[2], "` is missing.")
    abort(err_msg)
  }
}

