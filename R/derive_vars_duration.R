#' Derive Duration
#'
#' Derives duration between two dates, specified by the variables present in
#' input dataset e.g., duration of adverse events, relative day, age, ...
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @param new_var Name of variable to create
#'
#' @param new_var_unit Name of the unit variable If the parameter is not
#'   specified, no variable for the unit is created.
#'
#' @param start_date The start date
#'
#'   A date or date-time variable is expected. This variable must be present in
#'   specified input dataset.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param end_date The end date
#'
#'   A date or date-time variable is expected. This variable must be present in
#'   specified input dataset.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
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
#'   Default: `TRUE``
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @param add_one Add one input unit?
#'
#'   If the duration is non-negative, one input unit is added. I.e., the
#'   duration can not be zero.
#'
#'   Default: `TRUE` Permitted Values: `TRUE`, `FALSE`
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
#' @details The duration is derived as time from start to end date in the
#'   specified output unit. If the end date is before the start date, the duration
#'   is negative. The start and end date variable must be present in the specified
#'   input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the duration and unit variable added
#'
#' @keywords adam timing derivation
#'
#' @export
#'
#' @seealso [compute_duration()]
#'
#' @examples
#' library(lubridate)
#' library(tibble)
#'
#' # derive age in years
#' data <- tribble(
#'   ~BRTHDT, ~RANDDT,
#'   ymd("1984-09-06"), ymd("2020-02-24"),
#'   ymd("1985-01-01"), NA,
#'   NA, ymd("2021-03-10"),
#'   NA, NA
#' )
#'
#' derive_vars_duration(data,
#'   new_var = AAGE,
#'   new_var_unit = AAGEU,
#'   start_date = BRTHDT,
#'   end_date = RANDDT,
#'   out_unit = "years",
#'   add_one = FALSE,
#'   trunc_out = TRUE
#' )
#'
#' # derive adverse event duration in days
#' data <- tribble(
#'   ~ASTDT, ~AENDT,
#'   ymd("2021-03-05"), ymd("2021-03-02"),
#'   ymd("2019-09-18"), ymd("2019-09-18"),
#'   ymd("1985-01-01"), NA,
#'   NA, NA
#' )
#'
#' derive_vars_duration(data,
#'   new_var = ADURN,
#'   new_var_unit = ADURU,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   out_unit = "days"
#' )
#'
#' # derive adverse event duration in minutes
#' data <- tribble(
#'   ~ADTM,                         ~TRTSDTM,
#'   ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-09T05:00:00"),
#'   ymd_hms("2019-11-11T10:30:00"), ymd_hms("2019-11-11T11:30:00"),
#'   ymd("2019-11-11"),              ymd_hms("2019-11-11T04:00:00"),
#'   NA,                             ymd_hms("2019-11-11T12:34:56"),
#' )
#'
#' derive_vars_duration(data,
#'   new_var = ADURN,
#'   new_var_unit = ADURU,
#'   start_date = ADTM,
#'   end_date = TRTSDTM,
#'   in_unit = "minutes",
#'   out_unit = "minutes",
#'   add_one = FALSE
#' )
derive_vars_duration <- function(dataset,
                                 new_var,
                                 new_var_unit = NULL,
                                 start_date,
                                 end_date,
                                 in_unit = "days",
                                 out_unit = "days",
                                 floor_in = TRUE,
                                 add_one = TRUE,
                                 trunc_out = FALSE) {
  new_var <- assert_symbol(enquo(new_var))
  new_var_unit <- assert_symbol(enquo(new_var_unit), optional = TRUE)
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, required_vars = vars(!!start_date, !!end_date))
  assert_character_scalar(in_unit, values = valid_time_units())
  assert_character_scalar(out_unit, values = c(valid_time_units(), "weeks"))
  assert_logical_scalar(floor_in)
  assert_logical_scalar(add_one)
  assert_logical_scalar(trunc_out)

  warn_if_vars_exist(
    dataset,
    c(
      deparse(substitute(new_var)),
      deparse(substitute(new_var_unit))
    )
  )

  dataset <- dataset %>%
    mutate(
      !!new_var := compute_duration(
        !!start_date,
        !!end_date,
        in_unit = in_unit,
        out_unit = out_unit,
        floor_in = floor_in,
        add_one = add_one,
        trunc_out = trunc_out
      )
    )

  if (!quo_is_null(new_var_unit)) {
    dataset <- dataset %>%
      mutate(!!new_var_unit := if_else(is.na(!!new_var), NA_character_, toupper(out_unit)))
  }

  dataset
}
