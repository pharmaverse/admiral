#' Derive Duration
#'
#' Derives duration between two dates, e.g., duration of adverse events, relative
#' day, age, ...
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
#' @param end_date The end date
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
#' @param trunc_out Return integer part
#'
#'   The fractional part of the duration (in output unit) is removed, i.e., the
#'   integer part is returned.
#'
#'   Default: ``FALSE``
#'
#'   Permitted Values: ``TRUE``, ``FALSE``
#'
#' @details The duration is derived as time from start to end date in the
#'   specied output unit. If the end date is before the start date, the duration
#'   is negative.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the duration and unit variable added
#'
#' @family {general functions}
#'
#' @export
#'
#' @seealso [compute_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   ymd('1984-09-06'), ymd('2020-02-24'))
#'
#' derive_duration(data,
#'                 new_col = AAGE,
#'                 unit_col = AAGEU,
#'                 start_date = BRTHDT,
#'                 end_date = RANDDT,
#'                 out_unit = 'years',
#'                 add_one = FALSE,
#'                 trunc_out = TRUE)
#'


derive_duration <- function(dataset,
                            new_col,
                            unit_col,
                            start_date,
                            end_date,
                            in_unit = 'days',
                            out_unit = 'days',
                            floor_in = TRUE,
                            add_one = TRUE,
                            trunc_out = FALSE) {
  dataset <-
    dataset %>% mutate(!!enquo(new_col) := compute_duration(!!enquo(start_date),
                                                           !!enquo(end_date),
                                                           in_unit = in_unit,
                                                           out_unit = out_unit,
                                                           floor_in = floor_in,
                                                           add_one = add_one,
                                                           trunc_out = trunc_out))
  if(!missing(unit_col)){
    dataset <- dataset %>% mutate(!!enquo(unit_col) := toupper(out_unit))
  }
  return(dataset)
}