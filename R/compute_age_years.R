#' Compute Age in Years
#'
#' Converts a set of age values from the specified time unit to years.
#'
#' @param age The ages to convert.
#'
#'   A numeric vector is expected.
#'
#' @param age_unit Age unit.
#'
#'   Either a string containing the time unit of all ages in `age` or a character
#'   vector containing the time units of each age in `age` is expected. Note that
#'   permitted values are cases insensitive (e.g. `"YEARS"` is treated the same
#'   as `"years"` and `"Years"`).
#'
#'   Permitted Values: `"years"`, `"months"`, `"weeks"`, `"days"`, `"hours"`, `"minutes"`,
#'   `"seconds"`.
#'
#' @details Returns a numeric vector of ages in years as doubles. Note, underlying
#' computations assume an equal number of days in each year (365.25).
#'
#' @return The ages contained in `age` converted to years.
#'
#' @keywords com_date_time
#'
#' @family com_date_time
#'
#' @export
#'
#' @examples
#' compute_age_years(
#'   age = c(240, 360, 480),
#'   age_unit = "MONTHS"
#' )
#'
#' compute_age_years(
#'   age = c(10, 520, 3650),
#'   age_unit = c("YEARS", "WEEKS", "DAYS")
#' )
#'
compute_age_years <- function(age,
                              age_unit) {
  assert_numeric_vector(age)
  assert_character_vector(
    unique(tolower(age_unit)),
    values = c(
      NA, "years", "months", "weeks", "days",
      "hours", "minutes", "seconds"
    )
  )

  if (!(length(age_unit) %in% c(1, length(age)))) {
    abort(paste0(
      "`age_unit` must be a single string or a vector of the same length as",
      "`age`, but there are ", length(age), " values in `age` and ",
      length(age_unit), " values in `age_unit`."
    ))
  }

  age_years <- time_length(
    duration(age,
      units = tolower(age_unit)
    ),
    unit = "years"
  )

  age_years
}
