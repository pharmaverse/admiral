#' Derive Analysis End Relative Day
#'
#' Adds the analysis end relative day (`AENDY`) to the dataset, i.e. study day
#' of analysis end date
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @param start_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   The default is `TRTSDT`.
#'
#' @param end_date The end date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#'   The default is `AENDT`
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with `AENDY` column added
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~TRTSDT, ~AENDT,
#'   lubridate::ymd("2020-01-01"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_var_aendy(data)
derive_var_aendy <- function(dataset, start_date = TRTSDT, end_date = AENDT) {
  derive_duration(
    dataset,
    new_var = AENDY,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date)
  )
}
