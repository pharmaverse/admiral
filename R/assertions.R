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

######################################################
#assertions for date time derivations
######################################################
#check validity of the date imputation input
is_valid_date_entry<-function(arg){
  pattern<-"^([0-9]{2})-([0-9]{2})$"
  grepl(pattern,arg) | arg %in% c("FIRST", "MID", "LAST")
}
assertthat::on_failure(is_valid_date_entry) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid date entry.\n",
         "Date_imputation should be specified as 'dd-mm' (e.g. '01-01') or 'FIRST', 'MID', 'LAST' to get the first/mid/last day/month")
}
#check validity of the time imputation input
is_valid_time_entry<-function(arg){
  pattern<-"^([0-9]{2}):([0-9]{2}):([0-9]{2})$"
  grepl(pattern,arg) | arg %in% c("FIRST", "LAST")
}
assertthat::on_failure(is_valid_time_entry) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid time entry.\n",
         "time_imputation should be specified as 'hh:mm:ss' (e.g. '00:00:00') or 'FIRST','LAST' to get the first/last time of teh day")
}

#seconds, minutes, hours
is_valid_sec_min<-function(arg){
  arg %in% 0:59
}
assertthat::on_failure(is_valid_sec_min) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid min/sec.\n",
         "Values must be between between 0-59")
}

is_valid_hour<-function(arg){
  arg %in% 0:23
}
assertthat::on_failure(is_valid_hour) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid hour.\n",
         "Values must be between 0-23")
}

#day, month
is_valid_day<-function(arg){
  arg %in% 1:31 | arg %in% c("FIRST","MID", "LAST")
}
assertthat::on_failure(is_valid_day) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), "=", eval(call$arg, envir = env), " is not a valid day.\n",
         "Values must be between 1-31 or 'FISRT', 'MID', 'LAST'")
}


is_valid_month<-function(arg){
  arg %in% 1:12
}
assertthat::on_failure(is_valid_month) <- function(call, env) {
  paste0("Argument ", deparse(call$arg), " = ", eval(call$arg, envir = env), " is not a valid month.\n",
         "Values for month must be between 1-12. Please check the date_imputation input: it should be sepcified as 'dd-mm'")
}

is_valid_dtc<-function(arg){
  pattern0<-"^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2}):(\\d{2}).(\\d{3})$"
  pattern1<-"^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2}):(\\d{2})$"
  pattern2<-"^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2})$"
  pattern3<-"^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2})$"
  pattern4<-"^(\\d{4})-(\\d{2})-(\\d{2})$"
  pattern5<-"^(\\d{4})-(\\d{2})$"
  pattern6<-"^(\\d{4})$"
  pattern7<-"^(\\d{4})---(\\d{2})$"


  grepl(pattern0,arg)|
  grepl(pattern1,arg)|
  grepl(pattern2,arg)|
  grepl(pattern3,arg)|
  grepl(pattern4,arg)|
  grepl(pattern5,arg)|
  grepl(pattern6,arg)|
  grepl(pattern7,arg)
}

assert_is_valid_dtc <- function(dtc) {
  is_valid_dtc <- is_valid_dtc(dtc)

  if (!all(is_valid_dtc)) {
    incorrect_dtc <- dtc[is_valid_dtc==FALSE]
    incorrect_dtc_row <- rownames(as.data.frame(dtc))[is_valid_dtc==FALSE]
    tbl <- paste("Row: ",incorrect_dtc_row, ", --DTC: ",incorrect_dtc)
    msg<-"Dataset contains incorrect datetime format: --DTC may be incorrectly imputed on row(s)"
    warn(msg)
    warn(tbl)

    msg3 <- paste0(
      "The following representations are handled: \n",
      "2003-12-15T13:15:17.123\n",
      "2003-12-15T13:15:17\n",
      "2003-12-15T13:15\n",
      "2003-12-15T13\n",
      "2003-12-15\n",
      "2003-12\n",
      "2003\n",
      "2003---15\n\n",
      "The following representations are NOT handled: \n",
      "2003-12-15T-:15:18\n",
      "2003-12-15T13:-:19\n",
      "--12-15\n",
      "-----T07:15"
      )
    warn(msg3)

  }
}



