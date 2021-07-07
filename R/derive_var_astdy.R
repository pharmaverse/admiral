#' Derive Analysis Start Relative Day
#'
#' Adds the analysis start relative day (`ASTDY`) to the dataset, i.e., study
#' day of analysis start date.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `reference_date` and the `date` parameter are
#'   expected.
#'
#' @param reference_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   The default is `TRTSDT`.
#'
#' @param date The end date column for which the study day should be derived
#'
#'   A date or date-time object column is expected.
#'
#'   The default is `ASTDT`
#'
#' @author Stefan Bundfuss
#'
#' @details The study day is derived as number of days from the start date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   start date is 1.
#'
#' @return The input dataset with `ASTDY` column added
#'
#' @keywords adam timing derivation
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~TRTSDT, ~ASTDT,
#'   lubridate::ymd("2020-01-01"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_var_astdy(data)
derive_var_astdy <- function(dataset, reference_date = TRTSDT, date = ASTDT) {
  reference_date <- assert_symbol(enquo(reference_date))
  date <- assert_symbol(enquo(date))
  assert_data_frame(dataset, vars(!!reference_date, !!date))

  derive_duration(
    dataset,
    new_var = ASTDY,
    start_date = !!reference_date,
    end_date = !!date
  )
}
