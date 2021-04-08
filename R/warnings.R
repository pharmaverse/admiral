#' Warn If a Variable Already Exists
#'
#' Warn if a variable already exists inside a dataset
#'
#' @param dataset A `data.frame`
#' @param vars `character` vector of columns to check for in `dataset`
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @examples
#' data(dm)
#'
#' ## No warning as `AAGE` doesn't exist in `dm`
#' warn_if_vars_exist(dm, "AAGE")
#'
#' ## Issues a warning
#' warn_if_vars_exist(dm, "ARM")
warn_if_vars_exist <- function(dataset, vars) {
  existing_vars <- vars[vars %in% colnames(dataset)]
  if (length(existing_vars) == 1L) {
    msg <- paste("Variable", backquote(existing_vars), "already exists in the dataset")
    warn(msg)
  } else if (length(existing_vars) > 1L) {
    msg <- paste("Variables", enumerate(existing_vars), "already exist in the dataset")
    warn(msg)
  } else {
    invisible(NULL)
  }
}

is_valid_dtc <- function(arg) {
  pattern0 <- "^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2}):(\\d{2}).(\\d{3})$"
  pattern1 <- "^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2}):(\\d{2})$"
  pattern2 <- "^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2}):(\\d{2})$"
  pattern3 <- "^(\\d{4})-(\\d{2})-(\\d{2})T(\\d{2})$"
  pattern4 <- "^(\\d{4})-(\\d{2})-(\\d{2})$"
  pattern5 <- "^(\\d{4})-(\\d{2})$"
  pattern6 <- "^(\\d{4})$"
  pattern7 <- "^(\\d{4})---(\\d{2})$"


  grepl(pattern0, arg) |
    grepl(pattern1, arg) |
    grepl(pattern2, arg) |
    grepl(pattern3, arg) |
    grepl(pattern4, arg) |
    grepl(pattern5, arg) |
    grepl(pattern6, arg) |
    grepl(pattern7, arg) |
    arg == ""
}
#' Warn If a vector contains unknown datetime format
#'
#' Warn if the vector contains unknown datetime format such as
#' "2003-12-15T-:15:18", "2003-12-15T13:-:19","--12-15","-----T07:15"
#'
#' @param dtc a character vector containing the dates
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#'
#' ## No warning as `dtc` is a valid date format
#' warn_if_invalid_dtc(dtc = "2021-04-06")
#'
#' ## Issues a warning
#' warn_if_invalid_dtc(dtc = "2021-04-06T-:30:30")
warn_if_invalid_dtc <- function(dtc) {
  is_valid_dtc <- is_valid_dtc(dtc)

  if (!all(is_valid_dtc)) {
    incorrect_dtc <- dtc[is_valid_dtc == FALSE]
    incorrect_dtc_row <- rownames(as.data.frame(dtc))[is_valid_dtc == FALSE]
    tbl <- paste("Row", incorrect_dtc_row, ": --DTC = ", incorrect_dtc)
    msg <- "Dataset contains incorrect datetime format: --DTC may be incorrectly imputed on row(s)"
    warn(msg)
    warn(paste(capture.output(print(tbl)), collapse = "\n"))

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
