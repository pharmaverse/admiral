#' Derive Total Treatment Duration (Days)
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
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Default: `TRTSDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Default: `TRTEDT`
#'
#' @details The total treatment duration is derived as the number of days from
#'   start to end date plus one.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with `TRTDURD` added
#'
#' @family der_date_time
#'
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @seealso [derive_vars_duration()]
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#'
#' data <- tribble(
#'   ~TRTSDT, ~TRTEDT,
#'   ymd("2020-01-01"), ymd("2020-02-24")
#' )
#'
#' derive_var_trtdurd(data)
derive_var_trtdurd <- function(dataset,
                               start_date = TRTSDT,
                               end_date = TRTEDT) {
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, vars(!!start_date, !!end_date))

  derive_vars_duration(
    dataset,
    new_var = TRTDURD,
    start_date = !!start_date,
    end_date = !!end_date
  )
}
