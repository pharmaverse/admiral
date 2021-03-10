#' Derive Analysis End Relative Day
#'
#' Adds the analysis end relative day (``AENDY``) to the dataset, i.e., study day of analysis end date.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `startdate` and the `enddate` parameter are
#'   expected.
#'
#' @param startdate The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   The default is ``TRTSDT``.
#'
#' @param enddate The end date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#'   The default is ``AENDT``
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with ``AENDY`` column added
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~TRTSDT, ~AENDT,
#'   ymd('2020-01-01'), ymd('2020-02-24'))
#'
#' derive_var_aendy(data)

derive_var_aendy <- function(dataset, startdate = TRTSDT, enddate = AENDT){
  derive_duration(dataset,
                  newcol = AENDY,
                  startdate = !!enquo(startdate),
                  enddate = !!enquo(enddate))

}