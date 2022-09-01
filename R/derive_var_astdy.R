#' Derive Analysis Start Relative Day
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_dy()` instead.
#'
#' Adds the analysis start relative day (`ASTDY`) to the dataset, i.e., study
#' day of analysis start date.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `reference_date` and the `date` parameter are
#'   expected.
#'
#' @param reference_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character vector to a date object.
#'
#'   The default is `TRTSDT`.
#'
#' @param date The end date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character vector to a date object.
#'
#'   The default is `ASTDT`
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with `ASTDY` column added
#'
#' @keywords deprecated
#'
#' @export
#'
derive_var_astdy <- function(dataset, reference_date = TRTSDT, date = ASTDT) {
  deprecate_stop("0.7.0", "derive_var_astdy()", "derive_vars_dy()")
}
