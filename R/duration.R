#' Derive Duration
#'
#' Derive duration between two dates, e.g., duration of adverse events, relative
#' day, age, ...
#'
#' @param startdate The start date
#'
#'   A date or date-time object is expected.
#'
#' @param enddate The end date
#'
#'   A date or date-time object is expected.
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
#'   Permitted Values: 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
#'
#' @param floor_in Round down input dates?
#'
#'   The input dates are round down with respect to the input unit, e.g., if the
#'   input unit is 'days', the time of the input dates is ignored.
#'
#'   Default: ``TRUE```
#'
#'   Permitted Values: ``TRUE``, ``FALSE``
#'
#' @param add_one Add one input unit?
#'
#'   If the duration is non-negative, one input unit is added. I.e., the
#'   duration can not be zero.
#'
#'   Default: ``TRUE``
#'   Permitted Values: ``TRUE``, ``FALSE``
#'
#' @details The duration is derived as time from start to end date in the
#'   specied output unit. If the end date is before the start date, the duration
#'   is negative.
#'
#' @author Stefan Bundfuss
#'
#' @return The duration between the two date in the specified unit
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' # derive duration in days (integer), i.e., relative day
#' duration(ymd_hms('2020-12-06T15:00:00'), ymd_hms('2020-12-24T08:15:00'))
#'
#' # derive duration in days (float)
#' duration(ymd_hms('2020-12-06T15:00:00'), ymd_hms('2020-12-24T08:15:00'), floor_in = FALSE, add_one = FALSE)
#'
#' # derive age
#' duration(ymd('1984-09-06'), ymd('2020-12-24'), in_unit = 'years', out_unit = 'year', add_one = FALSE)
#'


duration <- function(startdate,
                     enddate,
                     in_unit = 'days',
                     out_unit = 'days',
                     floor_in = TRUE,
                     add_one = TRUE) {
    # Checks
    assert_that(arg_specified(startdate), arg_specified(enddate))
    assert_that(is_date(startdate), is_date(enddate))
    assert_that(is_timeunit(in_unit), is_timeunit(out_unit))

    # Derivation
    if (floor_in) {
      startdate <- floor_date(startdate, unit = in_unit)
      enddate <- floor_date(enddate, unit = in_unit)
    }

    duration <- time_length(startdate %--% enddate, unit = out_unit)
    if (add_one) {
      duration <- duration + time_length(if_else(duration >= 0,
                                                 lubridate::duration(1, unit = in_unit),
                                                 lubridate::duration(0, unit = in_unit)),
                                         unit = out_unit)
    }
    return(duration)
  }