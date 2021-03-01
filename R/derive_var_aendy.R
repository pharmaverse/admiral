#' Derive Analysis End Relative Day
#'
#' Adds the analysis end relative day (``AENDY``) to the dataset, i.e., study day of analysis end date.
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
#'   The default is ``AENDT``
#'
#' @inherit studyday details
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``AENDY`` column added
#'
#' @export
#'
#' @examples
#' derive_var_aendy(adae)

derive_var_aendy <- function(ds, refdate = TRTSDTM, date = AENDT){
  derive_studyday(ds, newcol = AENDY, refdate = !!enquo(refdate), date = !!enquo(date))
}