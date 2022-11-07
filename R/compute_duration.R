#' Compute Duration
#'
#' Compute duration between two dates, e.g., duration of an adverse event,
#' relative day, age, ...
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Refer to `convert_dtc_to_dt()` to obtain a vector of imputed dates.
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Refer to `convert_dtc_to_dt()` to obtain a vector of imputed dates.
#'
#' @param in_unit Input unit
#'
#'   See floor_in and add_one parameter for details.
#'
#'   Default: 'days'
#'
#'   Permitted Values: 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
#'
#' @param out_unit Output unit
#'
#'   The duration is derived in the specified unit
#'
#'   Default: 'days'
#'
#'   Permitted Values: 'years', 'months', 'weeks', 'days', 'hours', 'minutes',
#'   'seconds'
#'
#' @param floor_in Round down input dates?
#'
#'   The input dates are round down with respect to the input unit, e.g., if the
#'   input unit is 'days', the time of the input dates is ignored.
#'
#'   Default: `TRUE``
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @param add_one Add one input unit?
#'
#'   If the duration is non-negative, one input unit is added. i.e., the
#'   duration can not be zero.
#'
#'   Default: `TRUE`
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @param trunc_out Return integer part
#'
#'   The fractional part of the duration (in output unit) is removed, i.e., the
#'   integer part is returned.
#'
#'   Default: `FALSE`
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @details The output is a numeric vector providing the duration as time from
#' start to end date in the specified unit. If the end date is before the start
#' date, the duration is negative.
#'
#' @author Stefan Bundfuss
#'
#' @return The duration between the two date in the specified unit
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' # Derive duration in days (integer), i.e., relative day
#' compute_duration(
#'   start_date = ymd_hms("2020-12-06T15:00:00"),
#'   end_date = ymd_hms("2020-12-24T08:15:00")
#' )
#'
#' # Derive duration in days (float)
#' compute_duration(
#'   start_date = ymd_hms("2020-12-06T15:00:00"),
#'   end_date = ymd_hms("2020-12-24T08:15:00"),
#'   floor_in = FALSE,
#'   add_one = FALSE
#' )
#'
#' # Derive age in years
#' compute_duration(
#'   start_date = ymd("1984-09-06"),
#'   end_date = ymd("2020-02-24"),
#'   trunc_out = TRUE,
#'   out_unit = "years",
#'   add_one = FALSE
#' )
#'
#' # Derive duration in hours
#' compute_duration(
#'   start_date = ymd_hms("2020-12-06T9:00:00"),
#'   end_date = ymd_hms("2020-12-06T13:30:00"),
#'   out_unit = "hours",
#'   floor_in = FALSE,
#'   add_one = FALSE,
#' )
compute_duration <- function(start_date,
                             end_date,
                             in_unit = "days",
                             out_unit = "days",
                             floor_in = TRUE,
                             add_one = TRUE,
                             trunc_out = FALSE) {
  # Checks
  assert_date_vector(start_date)
  assert_date_vector(end_date)
  assert_character_scalar(in_unit, values = valid_time_units())
  assert_character_scalar(out_unit, values = c(valid_time_units(), "weeks"))
  assert_logical_scalar(floor_in)
  assert_logical_scalar(add_one)
  assert_logical_scalar(trunc_out)

  # Derivation
  if (floor_in) {
    # Remove information more precise than the input unit, e.g., if input unit
    # is days, the time part of the dates is removed. After updates in R `NA`
    # values have to be explicitly handled here because otherwise
    # `floor_date(as.Date(NA), unit = "days")` will return `"1970-01-01"` rather
    # than `NA` (#1486). See also here:
    # https://github.com/tidyverse/lubridate/issues/1069
    start_date_fun <- if (is.Date(start_date)) as.Date else as.POSIXct
    end_date_fun <- if (is.Date(end_date)) as.Date else as.POSIXct
    start_date <- if_else(
      is.na(start_date),
      start_date_fun(NA),
      floor_date(start_date, unit = in_unit)
    )
    end_date <- if_else(
      is.na(end_date),
      end_date_fun(NA),
      floor_date(end_date, unit = in_unit)
    )
  }

  # derive the duration in the output unit
  duration <- time_length(start_date %--% end_date, unit = out_unit)
  if (add_one) {
    # add one unit of the input unit (converted to the output unit), e.g., if
    # input unit is days and output unit is hours, 24 hours are added
    duration <- duration + time_length(
      if_else(
        duration >= 0,
        duration(1, unit = in_unit),
        duration(0, unit = in_unit)
      ),
      unit = out_unit
    )
  }
  if (trunc_out) {
    # remove fractional part
    duration <- trunc(duration)
  }
  duration
}
