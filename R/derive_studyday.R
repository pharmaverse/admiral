#' Derive Study Day
#'
#' Adds study day to dataset
#'
#' @param newcol Name of columnn to add
#'
#' @param refdate The reference date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#' @param date The date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#' @inherit studyday details
#'
#' @author Stefan Bundfuss
#'
#' @return The study day
#'
#' @export
#'
#' @examples
#' derive_studyday(adlb, ADY, refdate = TRTSDTM, date = ADT)


derive_studyday <- function(ds, newcol, refdate, date){
  # Checks

  # Derivation
  ds %>% mutate(!!enquo(newcol) := studyday(!!enquo(refdate), !!enquo(date)))
}
