#' Derive Analysis Age
#'
#' @description Derives analysis age (`AAGE`) and analysis age unit (`AAGEU`).
#'
#' **Note:** This is a wrapper function for the more generic `derive_vars_duration()`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("start_date", "end_date"))`
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character
#'   vector to a date object.
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character
#'   vector to a date object.
#'
#' @param age_unit Age unit
#'
#'   The age is derived in the specified unit
#'
#' @permitted The values are considered case-insensitive.
#'
#'   For years: `"year"`, `"years"`, `"yr"`, `"yrs"`, `"y"`
#'
#'   For months: `"month"`, `"months"`, `"mo"`, `"mos"`
#'
#'   For weeks: `"week"`, `"weeks"`, `"wk"`, `"wks"`, `"w"`
#'
#'   For days: `"day"`, `"days"`, `"d"`
#'
#'   For hours: `"hour"`, `"hours"`, `"hr"`, `"hrs"`, `"h"`
#'
#'   For minutes: `"minute"`, `"minutes"`, `"min"`, `"mins"`
#'
#'   For seconds: `"second"`, `"seconds"`, `"sec"`, `"secs"`, `"s"`
#'
#' @param type lubridate duration type
#'
#'   See below for details.
#'
#'   Default: `"interval"`
#'
#'   Permitted Values: `"duration"`, `"interval"`
#'
#' @details The duration is derived as time from start to end date in the
#'   specified output unit. If the end date is before the start date, the duration
#'   is negative. The start and end date variable must be present in the specified
#'   input dataset.
#'
#'   The [lubridate](https://lubridate.tidyverse.org/) package calculates two
#'   types of spans between two dates: duration and interval.
#'   While these calculations are largely the same, when the unit of the time period
#'   is month or year the result can be slightly different.
#'
#'   The difference arises from the ambiguity in the length of `"1 month"` or
#'   `"1 year"`.
#'   Months may have 31, 30, 28, or 29 days, and years are 365 days and 366 during leap years.
#'   Durations and intervals help solve the ambiguity in these measures.
#'
#'   The **interval** between `2000-02-01` and `2000-03-01` is `1` (i.e. one month).
#'   The **duration** between these two dates is `0.95`, which accounts for the fact
#'   that the year 2000 is a leap year, February has 29 days, and the average month
#'   length is `30.4375`, i.e. `29 / 30.4375 = 0.95`.
#'
#'   For additional details, review the
#'   [lubridate time span reference page](https://lubridate.tidyverse.org/reference/timespan.html).
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @family der_adsl
#' @keywords der_adsl
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
#'   ~BRTHDT, ~RANDDT,
#'   ymd("1984-09-06"), ymd("2020-02-24")
#' )
#'
#' derive_vars_aage(data)
derive_vars_aage <- function(dataset,
                             start_date = BRTHDT,
                             end_date = RANDDT,
                             age_unit = "YEARS",
                             type = "interval") {
  start_date <- assert_symbol(enexpr(start_date))
  end_date <- assert_symbol(enexpr(end_date))
  assert_data_frame(dataset, required_vars = expr_c(start_date, end_date))

  assert_character_scalar(age_unit, values = c(
    c("year", "years", "yr", "yrs", "y"),
    c("month", "months", "mo", "mos"),
    c("week", "weeks", "wk", "wks", "w"),
    c("day", "days", "d"),
    c("hour", "hours", "hr", "hrs", "h"),
    c("minute", "minutes", "min", "mins"),
    c("second", "seconds", "sec", "secs", "s")
  ), case_sensitive = FALSE)

  derive_vars_duration(
    dataset,
    new_var = AAGE,
    new_var_unit = AAGEU,
    start_date = !!start_date,
    end_date = !!end_date,
    out_unit = age_unit,
    add_one = FALSE,
    trunc_out = TRUE,
    type = type
  )
}


#' Derive Age in Years
#'
#' Converts the given age variable (`age_var`) to the unit 'years' from the current
#' units given in the `age_var+U` variable or `age_unit` argument and stores
#' in a new variable (`new_var`).
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("age_var"))`
#'
#' @param age_var Age variable.
#'
#'   A numeric object is expected.
#'
#' @param age_unit Age unit.
#'
#'   The `age_unit` argument is only expected when there is NOT a variable `age_var+U`
#'   in `dataset`. This gives the unit of the `age_var` variable and is used to convert
#'   AGE to 'years' so that grouping can occur.
#'
#'
#' @permitted 'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'
#'
#' @param new_var New age variable to be created in years. The returned values are
#'   doubles and NOT integers.
#''
#' @details This function is used to convert an age variable into the unit 'years'
#'   which can then be used to create age groups. The resulting column contains the
#'   equivalent years as a double. Note, underlying computations assume an equal number
#'   of days in each year (365.25).
#'
#' @return The input dataset (`dataset`) with `new_var` variable added in years.
#'
#' @family der_adsl
#' @keywords der_adsl
#'
#' @export
#'
#' @seealso [derive_vars_duration()]
#'
#' @examples
#' library(tibble)
#'
#' # Derive age with age units specified
#' data <- tribble(
#'   ~AGE, ~AGEU,
#'   27, "days",
#'   24, "months",
#'   3, "years",
#'   4, "weeks",
#'   1, "years"
#' )
#'
#' derive_var_age_years(data, AGE, new_var = AAGE)
#'
#' # Derive age without age units variable specified
#' data <- tribble(
#'   ~AGE,
#'   12,
#'   24,
#'   36,
#'   48
#' )
#' derive_var_age_years(data, AGE, age_unit = "months", new_var = AAGE)
derive_var_age_years <- function(dataset, age_var, age_unit = NULL, new_var) {
  age_variable <- assert_symbol(enexpr(age_var))
  assert_data_frame(dataset, required_vars = expr_c(age_variable))

  age_var <- pull(dataset, !!age_variable)
  assert_numeric_vector(age_var)

  age_var <- age_variable
  unit_var <- paste0(age_var, "U")

  age_unit <- assert_character_scalar(
    age_unit,
    values = c(
      "years", "months", "weeks", "days",
      "hours", "minutes", "seconds"
    ),
    case_sensitive = FALSE,
    optional = TRUE
  )

  new_var <- assert_symbol(enexpr(new_var))
  warn_if_vars_exist(dataset, as_name(new_var))

  if (!unit_var %in% colnames(dataset)) {
    if (is.null(age_unit)) {
      cli_abort(paste(
        "There is no unit variable ({.var {unit_var}}) associated with {.var {age_var}}.",
        "Please specify a value for {.arg age_unit}."
      ))
    } else {
      ds <- dataset %>%
        mutate(!!new_var := compute_age_years(!!age_var, age_unit))
    }
  } else {
    unit <- unique(tolower(pull(dataset, !!sym(unit_var))))
    assert_character_vector(
      unit,
      values = c(
        NA, "years", "months", "weeks", "days",
        "hours", "minutes", "seconds"
      )
    )

    if (!is.null(age_unit)) {
      if (length(unit) > 1) {
        msg <- paste(
          "The unit variable {.var {unit_var}} is associated with {.var {age_var}}",
          "and contains multiple values but the argument {.arg age_unit}",
          "has been specified with a single different value.",
          "The {.arg age_unit} argument is ignored and the conversion will be based",
          "on {.var {unit_var}}."
        )
        cli_warn(msg)
      } else if (unit != age_unit) {
        msg <- paste(
          "The unit variable {.var {unit_var}} is associated with {.var {age_var}}",
          "but the argument {.arg age_unit} has been specified with a different value.",
          "The {.arg age_unit} argument is ignored and the conversion will be based",
          "on {.var {unit_var}}."
        )
        cli_warn(msg)
      }
    }

    ds <- dataset %>%
      mutate(!!new_var := compute_age_years(!!age_var, !!sym(unit_var)))
  }
  return(ds)
}
