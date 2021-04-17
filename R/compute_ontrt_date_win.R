#' Compute date with window added
#'
#' @param ref_date Date variable for calculation to be done on; Required
#'
#' @param ref_window Amount of time to add to `ref_data`.  The unit is
#' specified in `ref_window_units`; Required
#'
#' @param ref_window_units A time unit for `ref_window`.
#'   Valid values: days, weeks, years; required
#'
#' @details
#' Returns a date with a specified period added
#'
#' @author Alice Ehmann
#'
#' @return A date with window added
#'
#' @export
#'
#' @examples
#' ~ref_date, ~ref_window, ~ref_window_units, ~result
#' 2020-01-01, 1, days, 2020-01-02
#' 2020-01-01, 2, years, 2021-01-01
#' 2020-01-01, 1, weeks, 2020-01-08
#'
#' result = compute_ontrt_date_win(ref_date, 1, "days")
#' result = compute_ontrt_date_win(ref_date, 2, "years")
#' result = compute_ontrt_date_win(ref_date, 1, "weeks")
#'
#'
compute_ontrt_date_win <- function(ref_date,
                            ref_window,
                            ref_window_units) {

  # Checks
  assert_that(is_date(ref_date))
  assert_that(is.numeric(ref_window))
  arg_match(ref_window_units, c("days", "weeks", "years"))

  # Derive Date with window
  if (ref_window_units == "days") {
    ref_window_date <- ref_date + days(x = ref_window)
  } else if (ref_window_units == "weeks") {
    ref_window_date <- ref_date + weeks(x = ref_window)
  } else if (ref_window_units == "years") {
    ref_window_date <- ref_date + years(x = ref_window)
  }

  ref_window_date

}
