#' Derive Analysis Study Day
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_dy()` instead.
#'
#' Adds the analysis study day (`ADY`) to the dataset, i.e., study day of analysis date.
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
#'   The default is `ADT`
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is non-negative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with `ADY` column added
#'
#' @keywords deprecated
#'
#' @export
#'
derive_var_ady <- function(dataset, reference_date = TRTSDT, date = ADT) {
  deprecate_stop("0.8.0", "derive_var_ady()", "derive_vars_dy()")
}
