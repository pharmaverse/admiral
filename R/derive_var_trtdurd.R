#' Derive Total Treatment Duration (days)
#'
#' Derives total treatment duration (days) (`TRTDURD`)
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `startdate` and the `enddate` parameter are
#'   expected.
#'
#' @param startdate The start date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `TRTSDT`
#'
#' @param enddate The end date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `TRTEDT`
#'
#' @details The total treatment duration is derived as the number of days from start to
#'   end date plus one.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with `TRTDURD` added
#'
#' @export
#'
#' @seealso [derive_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~TRTSDT, ~TRTEDT,
#'   ymd('2020-01-01'), ymd('2020-02-24'))
#'
#' derive_var_trtdurd(data)
#'

derive_var_trtdurd <- function(dataset,
                               startdate = TRTSDT,
                               enddate = TRTEDT){
  derive_duration(dataset,
                  newcol = TRTDURD,
                  startdate = !!enquo(startdate),
                  enddate = !!enquo(enddate))
}