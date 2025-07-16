#' Derive Duration
#'
#' Derives duration between two dates, specified by the variables present in
#' input dataset e.g., duration of adverse events, relative day, age, ...
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("start_date", "end_date"))`
#'
#' @param new_var Name of variable to create
#'
#' @param new_var_unit Name of the unit variable If the parameter is not
#'   specified, no variable for the unit is created.
#'
#' @inheritParams compute_duration
#'
#' @details The duration is derived as time from start to end date in the
#'   specified output unit. If the end date is before the start date, the duration
#'   is negative. The start and end date variable must be present in the specified
#'   input dataset.
#'
#' The [lubridate](https://lubridate.tidyverse.org/) package calculates two
#' types of spans between two dates: duration and interval.
#' While these calculations are largely the same, when the unit of the time period
#' is month or year the result can be slightly different.
#'
#' The difference arises from the ambiguity in the length of `"1 month"` or
#' `"1 year"`.
#' Months may have 31, 30, 28, or 29 days, and years are 365 days and 366 during leap years.
#' Durations and intervals help solve the ambiguity in these measures.
#'
#' The **interval** between `2000-02-01` and `2000-03-01` is `1` (i.e. one month).
#' The **duration** between these two dates is `0.95`, which accounts for the fact
#' that the year 2000 is a leap year, February has 29 days, and the average month
#' length is `30.4375`, i.e. `29 / 30.4375 = 0.95`.
#'
#' For additional details, review the
#' [lubridate time span reference page](https://lubridate.tidyverse.org/reference/timespan.html).
#'
#'
#' @return The input dataset with the duration and unit variable added
#'
#' @family der_date_time
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @seealso [compute_duration()]
#'
#' @examples
#' library(lubridate)
#' library(tibble)
#'
#' # Derive age in years
#' data <- tribble(
#'   ~USUBJID, ~BRTHDT, ~RANDDT,
#'   "P01", ymd("1984-09-06"), ymd("2020-02-24"),
#'   "P02", ymd("1985-01-01"), NA,
#'   "P03", NA, ymd("2021-03-10"),
#'   "P04", NA, NA
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
#' # Derive adverse event duration in days
#' data <- tribble(
#'   ~USUBJID, ~ASTDT, ~AENDT,
#'   "P01", ymd("2021-03-05"), ymd("2021-03-02"),
#'   "P02", ymd("2019-09-18"), ymd("2019-09-18"),
#'   "P03", ymd("1985-01-01"), NA,
#'   "P04", NA, NA
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
#' # Derive adverse event duration in minutes
#' data <- tribble(
#'   ~USUBJID, ~ADTM, ~TRTSDTM,
#'   "P01", ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-09T05:00:00"),
#'   "P02", ymd_hms("2019-11-11T10:30:00"), ymd_hms("2019-11-11T11:30:00"),
#'   "P03", ymd_hms("2019-11-11T00:00:00"), ymd_hms("2019-11-11T04:00:00"),
#'   "P04", NA, ymd_hms("2019-11-11T12:34:56"),
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
#'
#' # Derive adverse event start time since last dose in hours
#' data <- tribble(
#'   ~USUBJID, ~ASTDTM, ~LDOSEDTM,
#'   "P01", ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-08T10:05:00"),
#'   "P02", ymd_hms("2019-11-11T23:59:59"), ymd_hms("2019-10-11T11:37:00"),
#'   "P03", ymd_hms("2019-11-11T00:00:00"), ymd_hms("2019-11-10T23:59:59"),
#'   "P04", ymd_hms("2019-11-11T12:34:56"), NA,
#'   "P05", NA, ymd_hms("2019-09-28T12:34:56")
#' )
#' derive_vars_duration(
#'   data,
#'   new_var = LDRELTM,
#'   new_var_unit = LDRELTMU,
#'   start_date = LDOSEDTM,
#'   end_date = ASTDTM,
#'   in_unit = "hours",
#'   out_unit = "hours",
#'   add_one = FALSE
#' )
derive_vars_duration <- function(dataset,
                                 new_var,
                                 new_var_unit = NULL,
                                 start_date,
                                 end_date,
                                 in_unit = "days",
                                 out_unit = "DAYS",
                                 floor_in = TRUE,
                                 add_one = TRUE,
                                 trunc_out = FALSE,
                                 type = "duration") {
  new_var <- assert_symbol(enexpr(new_var))
  new_var_unit <- assert_symbol(enexpr(new_var_unit), optional = TRUE)
  start_date <- assert_symbol(enexpr(start_date))
  end_date <- assert_symbol(enexpr(end_date))
  assert_data_frame(dataset, required_vars = exprs(!!start_date, !!end_date))

  in_unit <- get_unified_time_unit(in_unit)
  original_out_unit <- out_unit
  out_unit <- get_unified_time_unit(out_unit)

  assert_character_scalar(in_unit, values = c(
    c("year", "years", "yr", "yrs", "y"),
    c("month", "months", "mo", "mos"),
    c("day", "days", "d"),
    c("hour", "hours", "hr", "hrs", "h"),
    c("minute", "minutes", "min", "mins"),
    c("second", "seconds", "sec", "secs", "s")
  ))
  assert_character_scalar(out_unit, values = c(
    c("year", "years", "yr", "yrs", "y"),
    c("month", "months", "mo", "mos"),
    c("week", "weeks", "wk", "wks", "w"),
    c("day", "days", "d"),
    c("hour", "hours", "hr", "hrs", "h"),
    c("minute", "minutes", "min", "mins"),
    c("second", "seconds", "sec", "secs", "s")
  ))
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
        trunc_out = trunc_out,
        type = type
      )
    )

  if (!is.null(new_var_unit)) {
    dataset <- dataset %>%
      mutate(!!new_var_unit := if_else(is.na(!!new_var), NA_character_, original_out_unit))
  }

  dataset
}
