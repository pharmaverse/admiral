#' Derive Study Day
#'
#' Derive study day
#'
#' @param refdate The reference date, e.g., date of first treatment
#'
#'   A date or date-time object is expected.
#'
#' @param date The date for which the study day should be derived
#'
#'   A date or date-time object is expected.
#'
#' @details The study day is derived as number of days from the reference date
#'   to the date. If it is nonnegative, one is added. I.e., the study day of the
#'   reference date is 1.
#'
#' @author Stefan Bundfuss
#'
#' @return The study day
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' studyday(ymd('20201206'), ymd('20201224'))


studyday <- function(refdate, date){
  # Checks
  assert_that(arg_specified(refdate), arg_specified(date))
  assert_that(is_date(refdate), is_date(date))

  # Derivation
  studyday <- as_date(refdate) %--% as_date(date) / ddays()
  studyday <- studyday + ifelse(studyday>=0, 1, 0)

  return(studyday)
}
