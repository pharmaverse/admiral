#' Derive Total Treatment Duration (Days)
#'
#' @description Derives total treatment duration (days) (`TRTDURD`).
#'
#' **Note:** This is a wrapper function for the more generic `derive_vars_duration()`.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("start_date", "end_date"))`
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'
#' @details The total treatment duration is derived as the number of days from
#'   start to end date plus one.
#'
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
  start_date <- assert_symbol(enexpr(start_date))
  end_date <- assert_symbol(enexpr(end_date))
  assert_data_frame(dataset, exprs(!!start_date, !!end_date))

  derive_vars_duration(
    dataset,
    new_var = TRTDURD,
    start_date = !!start_date,
    end_date = !!end_date
  )
}
