#' Derive Analysis Study Day
#'
#' Adds the analysis study day (``ADY``) to the dataset, i.e., study day of analysis date.
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
#'   The default is ``ADT``
#'
#' @inherit studyday details
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``ADY`` column added
#'
#' @export
#'
#' @examples
#' derive_var_ady(adlb)

derive_var_ady <- function(ds, refdate = TRTSDTM, date = ADT){
  derive_studyday(ds, newcol = ADY, refdate = !!enquo(refdate), date = !!enquo(date))
}