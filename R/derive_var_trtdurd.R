#' Derive Total Treatment Duration (days)
#'
#' Derives total treatment duration (days) (`TRTDURD`)
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `TRTSDT`
#'
#' @param end_date The end date
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
                               start_date = TRTSDT,
                               end_date = TRTEDT){
  derive_duration(dataset,
                  newcol = TRTDURD,
                  start_date = !!enquo(start_date),
                  end_date = !!enquo(end_date))
}