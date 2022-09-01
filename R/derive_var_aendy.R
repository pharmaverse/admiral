#' Derive Analysis End Relative Day
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_dy()` instead.
#'
#' Adds the analysis end relative day (`AENDY`) to the dataset, i.e. study day
#' of analysis end date
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
#'   The default is `AENDT`
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with `AENDY` column added
#'
#' @keywords deprecated
#' @family deprecated
#' @export
#'
derive_var_aendy <- function(dataset, reference_date = TRTSDT, date = AENDT) {
  deprecate_stop("0.8.0", "derive_var_aendy()", "derive_vars_dy()")
}
