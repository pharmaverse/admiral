#' Derive Analysis Age
#'
#' Derives analysis age (`AAGE`) and analysis age unit (`AAGEU`)
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
#'   Default: `BRTHDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `RANDDT`
#'
#' @param unit Unit
#'
#'   The age is derived in the specified unit
#'
#'   Default: 'years'
#'
#'   Permitted Values: 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
#'
#' @details The age is derived as the integer part of the duration from start to
#'   end date in the specified unit.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @export
#'
#' @seealso [derive_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_aage(data)
derive_aage <- function(dataset,
                        start_date = BRTHDT,
                        end_date = RANDDT,
                        unit = "years") {
  derive_duration(
    dataset,
    new_var = AAGE,
    new_var_unit = AAGEU,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date),
    out_unit = unit,
    add_one = FALSE,
    trunc_out = TRUE
  )
}
