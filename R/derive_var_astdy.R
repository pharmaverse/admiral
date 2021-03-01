#' Derive Analysis Start Relative Day
#'
#' Adds the analysis start relative day (``ASTDY``) to the dataset, i.e., study day of analysis start date.
#'
#' @param refdate The reference date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   The default is ``TRTSDTM``.
#'
#' @param date The date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#'   The default is ``ASTDT``
#'
#' @inherit studyday details
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``ASTDY`` column added
#'
#' @export
#'
#' @examples
#' derive_var_astdy(adae)

derive_var_astdy <- function(ds, refdate = TRTSDTM, date = ASTDT){
  derive_studyday(ds, newcol = ASTDY, refdate = !!enquo(refdate), date = !!enquo(date))
}