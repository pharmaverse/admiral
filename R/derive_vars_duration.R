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
#' @param new_var Name of variable to create
#'
#' @param new_var_unit Name of the unit variable
#'   If the parameter is not specified, no variable for the unit is created.
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
#'   Default: `TRUE``
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @param add_one Add one input unit?
#'
#'   If the duration is non-negative, one input unit is added. I.e., the
#'   duration can not be zero.
#'
#'   Default: `TRUE`
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
#' @details The duration is derived as time from start to end date in the
#'   specified output unit. If the end date is before the start date, the duration
#'   is negative.
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
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
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
  assert_character_scalar(out_unit, values = valid_time_units())
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
    dataset <- dataset %>% mutate(!!new_var_unit := toupper(out_unit))
  }

  dataset
}

#' Derive Duration
#'
#' `derive_duration()` was renamed to `derive_vars_duration()` to create a
#' more consistent API.
#'
#' @keywords internal
#'
#' @export
derive_duration <- function(dataset,
                            new_var,
                            new_var_unit = NULL,
                            start_date,
                            end_date,
                            in_unit = "days",
                            out_unit = "days",
                            floor_in = TRUE,
                            add_one = TRUE,
                            trunc_out = FALSE) {
  deprecate_warn("0.3.0", "derive_duration()", "derive_vars_duration()")
  derive_vars_duration(
    dataset,
    new_var = !!enquo(new_var),
    new_var_unit = !!enquo(new_var_unit),
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date),
    in_unit = in_unit,
    out_unit = out_unit,
    floor_in = floor_in,
    add_one = add_one,
    trunc_out = trunc_out
  )
}
