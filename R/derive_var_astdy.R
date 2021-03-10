#' Derive Analysis Start Relative Day
#'
#' Adds the analysis start relative day (``ASTDY``) to the dataset, i.e., study
#' day of analysis start date.
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
#'   The default is ``ASTDT``
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with ``ASTDY`` column added
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~TRTSDT, ~ASTDT,
#'   ymd('2020-01-01'), ymd('2020-02-24'))
#'
#' derive_var_astdy(data)

derive_var_astdy <- function(dataset, startdate = TRTSDT, enddate = ASTDT){
  derive_duration(dataset,
                  newcol = ASTDY,
                  startdate = !!enquo(startdate),
                  enddate = !!enquo(enddate))
}