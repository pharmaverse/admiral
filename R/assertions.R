#' Argument Specified
#'
#' Checks if an argument was specified
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return ``TRUE`` if the argument was specified, ``FALSE`` otherwise
#'
#' @export
#'
#' @examples
#' assert_that(arg_specified(refdate), arg_specified(date))


arg_specified <- function(arg){
  return(!missing(arg))
}

assertthat::on_failure(arg_specified) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " is not specified.")
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
  return(is.instant(arg))
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
  return(arg %in% c('years', 'months', 'days', 'hours', 'minutes', 'seconds'))
}

assertthat::on_failure(is_timeunit) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid time unit.",
         " Valid time units are 'years', 'months', 'days', 'hours', 'minutes', and 'seconds'.")
}
