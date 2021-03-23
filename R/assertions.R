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
      err_msg <- paste0("Required variables ", enumerate(missing_vars), " are missing.")
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
#' assert_that(is_date(refdate), is_date(date))
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
#' assert_that(is_timeunit(unit))
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
