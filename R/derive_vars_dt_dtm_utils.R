#' Create a `dtm_level` object
#'
#' @param level Datetime level
#'
#' @permitted `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"h"` (hour), `"m"` (minute), `"s"` (second, lowest level), `"n"`
#'   (none)
#'
#' @returns A `dtm_level` object
#'
#' @details A `dtm_level` object is an ordered factor, i.e., two objects can be
#'   compared.
#'
#' @examples
#' # Create a dtm_level object with level "D" (day)
#' level_day <- dtm_level("D")
#' print(level_day)
#'
#' # Create a dtm_level object with level "h" (hour)
#' level_hour <- dtm_level("h")
#' print(level_hour)
#'
#' # Compare two dtm_level objects
#' level_day > level_hour # TRUE, because "D" is larger than "h".
#' @family utils_impute
#' @keywords internal
dtm_level <- function(level) {
  possible_values <- c("n", "s", "m", "h", "D", "M", "Y")
  assert_character_scalar(level, values = possible_values)

  out <-
    factor(
      level,
      levels = possible_values,
      ordered = TRUE
    )
  class(out) <- c("dtm_level", class(out))
  out
}

#' Create a `dt_level` object
#'
#' @param level Date level
#'
#' @permitted `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"n"` (none, lowest level)
#'
#' @returns A `dt_level` object
#'
#' @details A `dt_level` object is an ordered factor, i.e., two objects can be
#'   compared.
#'
#' @examples
#' # Create a dt_level object with level "D" (day)
#' level_day <- admiral:::dt_level("D")
#' print(level_day)
#'
#' # Create a dt_level object with level "Y" (year)
#' level_year <- admiral:::dt_level("Y")
#' print(level_year)
#'
#' # Compare two dt_level objects
#' level_day > level_year # TRUE, because "Y" is larger than "D".
#'
#' @family utils_impute
#' @keywords internal
dt_level <- function(level) {
  possible_values <- c("n", "D", "M", "Y")
  assert_character_scalar(level, values = possible_values)

  out <-
    factor(
      level,
      levels = possible_values,
      ordered = TRUE
    )
  class(out) <- c("dt_level", class(out))
  out
}

#' Get Date Imputation Targets
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for the 15th
#'   of June,
#'   - or as a keyword: `"first"`, `"mid"`, `"last"` to impute to the first/mid/last
#'   day/month.
#'
#' @param month Month component of the partial date
#'
#' @returns A list of character vectors. The elements of the list are named
#'   "year", "month", "day".
#'
#' @details
#'
#'  - For `date_imputation = "first"` `"0000"`, `"01"`, `"01"` are returned.
#'  - For `date_imputation = "mid"` `"xxxx"`, `"06"`, `"30"` if `month` is `NA`
#'  and `"15"` otherwise are returned.
#'  - For `date_imputation = "last"` `"9999"`, `"12"`, `"31"` are returned.
#'  - For `date_imputation = "<mm>-<dd>"` `"xxxx"`, `"<mm>"`, `"<dd>"` are returned.
#'
#'  `"xxxx"` indicates that the component is undefined. If an undefined
#'  component occurs in the imputed DTC value, the imputed DTC value is set to
#'  `NA_character_` in the imputation functions.
#'
#' @examples
#' # Get imputation target for "first"
#' target_first <- get_imputation_target_date("first", month = NA)
#' print(target_first)
#'
#' # Get imputation target for "mid" with specified month
#' target_mid <- get_imputation_target_date("mid", month = "03")
#' print(target_mid)
#'
#' # Get imputation target for "mid" with NA month
#' target_mid_na <- get_imputation_target_date("mid", month = NA)
#' print(target_mid_na)
#'
#' # Get imputation target for "last"
#' target_last <- get_imputation_target_date("last", month = NA)
#' print(target_last)
#'
#' # Get imputation target for custom date imputation "06-15"
#' target_custom <- get_imputation_target_date("06-15", month = NA)
#' print(target_custom)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_imputation_target_date <- function(date_imputation,
                                       month) {
  target <- vector("list", 3)
  names(target) <- c("year", "month", "day")
  target <- case_when(
    date_imputation == "first" ~ list(
      year = "0000",
      month = "01",
      day = "01"
    ),
    date_imputation == "mid" ~ list(
      year = "xxxx",
      month = "06",
      day = if_else(is.na(month), "30", "15")
    ),
    date_imputation == "last" ~ list(
      year = "9999",
      month = "12",
      day = "28"
    ),
    TRUE ~ list(
      year = "xxxx",
      month = str_sub(date_imputation, 1, 2),
      day = str_sub(date_imputation, 4, 5)
    )
  )

  target
}

#' Get Time Imputation Targets
#'
#' @param time_imputation The value to impute the time when a timepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as `"hh:mm:ss"`: e.g. `"00:00:00"`
#'   for the start of the day,
#'   - or as a keyword: `"first"`,`"last"` to impute to the start/end of a day.
#'
#' @returns A list of character vectors. The elements of the list are named
#'   "hour", "minute", "second".
#'
#' @details
#'
#'  - For `time_imputation = "first"` `"00"`, `"00"`, `"00"` are returned.
#'  - For `time_imputation = "last"` `"23"`, `"59"`, `"59"` are returned.
#'  - For `time_imputation = "<hh>:<mm>:<ss>"` `"<hh>"`, `"<mm>"`, `"<ss>"` are returned.
#'
#' @examples
#' # Get imputation target for "first" time
#' target_first_time <- get_imputation_target_time("first")
#' print(target_first_time)
#'
#' # Get imputation target for "last" time
#' target_last_time <- get_imputation_target_time("last")
#' print(target_last_time)
#'
#' # Get imputation target for custom time imputation "12-34-56"
#' target_custom_time <- get_imputation_target_time("12-34-56")
#' print(target_custom_time)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso  [impute_dtc_dtm()]
get_imputation_target_time <- function(time_imputation) {
  target <- vector("list", 3)
  names(target) <- c("hour", "minute", "second")
  target <- case_when(
    time_imputation == "first" ~ list(hour = "00", minute = "00", second = "00"),
    time_imputation == "last" ~ list(hour = "23", minute = "59", second = "59"),
    TRUE ~ list(
      hour = str_sub(time_imputation, 1, 2),
      minute = str_sub(time_imputation, 4, 5),
      second = str_sub(time_imputation, 7, -1)
    )
  )

  target
}

#' Convert a Date into a Datetime Object
#'
#' @description Convert a date (datetime, date, or date character) into a Date
#' vector (usually `'--DTM'`).
#'
#' **Note:** This is a wrapper function for the function `convert_dtc_to_dtm()`.
#'
#' @param dt The date to convert.
#'
#'   A date or character date is expected in a format like `yyyy-mm-ddThh:mm:ss`.
#'
#' @inheritParams convert_dtc_to_dtm
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @returns A datetime object
#'
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' convert_date_to_dtm("2019-07-18T15:25:00")
#' convert_date_to_dtm(Sys.time())
#' convert_date_to_dtm(as.Date("2019-07-18"), time_imputation = "23:59:59")
#' convert_date_to_dtm("2019-07-18", time_imputation = "23:59:59")
#' convert_date_to_dtm("2019-07-18")
convert_date_to_dtm <- function(dt,
                                highest_imputation = "h",
                                date_imputation = "first",
                                time_imputation = "first",
                                min_dates = NULL,
                                max_dates = NULL,
                                preserve = FALSE) {
  if (is.POSIXct(dt)) {
    return(dt)
  } else {
    if (is.instant(dt)) {
      dt <- format(dt, "%Y-%m-%d")
    }

    # convert dtc to dtm
    dt %>%
      convert_dtc_to_dtm(
        highest_imputation = highest_imputation,
        date_imputation = date_imputation,
        time_imputation = time_imputation,
        min_dates = min_dates,
        max_dates = max_dates,
        preserve = preserve
      )
  }
}

#' Parse DTC variable and Determine Components
#'
#' @param dtc The `'--DTC'` date to parse
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
#'   valid value for any component.
#'
#' @returns A list of character vectors. The elements of the list are named
#'   "year", "month", "day", "hour", "minute", and "second". Missing components
#'   are set to `NA_character_`.
#'
#' @details The function can be replaced by the parttime parser once it is
#'   available.
#'
#' @examples
#' # Get partial datetime components for a complete datetime string
#' dtc_complete <- get_partialdatetime("2020-03-15T12:34:56")
#' print(dtc_complete)
#'
#' # Get partial datetime components for a partial datetime string
#' dtc_partial <- get_partialdatetime("2020-03-15T12:34")
#' print(dtc_partial)
#'
#' # Get partial datetime components for a date-only string
#' dtc_date_only <- get_partialdatetime("2020-03-15")
#' print(dtc_date_only)
#'
#' # Get partial datetime components for an incomplete year string
#' dtc_year_partial <- get_partialdatetime("2020")
#' print(dtc_year_partial)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_partialdatetime <- function(dtc) {
  two <- "(\\d{2}|-?)"
  partialdate <- str_match(dtc, paste0(
    "(\\d{4}|-?)-?",
    two,
    "-?",
    two,
    "T?",
    two,
    ":?",
    two,
    ":?",
    "(\\d{2}(\\.\\d{1,5})?)?"
  ))
  partial <- vector("list", 6)
  components <- c("year", "month", "day", "hour", "minute", "second")
  names(partial) <- components
  for (i in seq_along(components)) {
    partial[[i]] <- partialdate[, i + 1]
    partial[[i]] <- if_else(partial[[i]] %in% c("-", ""), NA_character_, partial[[i]])
  }
  partial
}

#' Check input assertions for date and imputation parameters
#' Used in `derive_vars_dt()` and `derive_vars_dtm()`.
#'
#' Applies assertions on the `new_vars_prefix`, `max_dates`, `min_dates`,
#' `flag_imputation`, `highest_imputation`, and `date_imputation` arguments
#' to reduce cyclomatic complexity.
#'
#' @param new_vars_prefix Prefix for new variables.
#' @param max_dates Maximum dates of events (can be `NULL`).
#' @param min_dates Minimum dates of events (can be `NULL`).
#' @param flag_imputation_values Allowed values for the `flag_imputation` argument.
#' @param flag_imputation The value to impute.
#' @param highest_imputation Highest imputation level.
#' @param highest_imputation_values Allowed values for the `highest_imputation` argument.
#' @param date_imputation The value to impute the day/month
#' when a date part is missing (can be `NULL`).
#'
#' @keywords internal
#'
#' @returns `invisible(NULL)`
assert_dt_dtm_inputs <- function(new_vars_prefix, max_dates, min_dates,
                                 flag_imputation, flag_imputation_values,
                                 highest_imputation, highest_imputation_values,
                                 date_imputation = NULL) {
  assert_character_scalar(new_vars_prefix)

  assert_character_scalar(
    flag_imputation,
    values = flag_imputation_values,
    case_sensitive = FALSE
  )

  assert_highest_imputation(
    highest_imputation = highest_imputation,
    highest_imputation_values = highest_imputation_values,
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates
  )

  return(invisible(NULL))
}

#' Assert `date_imputation`
#'
#' Applies assertions on the `date_imputation` argument to reduce
#' cyclomatic complexity
#'
#' @param highest_imputation Highest imputation level
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#' @returns asserted `date_imputation`
#'
#' @keywords internal
#'
#' @returns `invisible(NULL)`
assert_date_imputation <- function(highest_imputation, date_imputation) {
  date_imputation <-
    assert_character_scalar(
      date_imputation,
      case_sensitive = FALSE
    )

  date_imputation <- tolower(date_imputation)
  if (highest_imputation == "D") {
    assert_character_scalar(date_imputation, values = c("first", "mid", "last"))
  }

  if (highest_imputation == "M") {
    is_mm_dd_format <- grepl("^(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$", date_imputation)
    is_one_of_keys <- date_imputation %in% c("first", "mid", "last")
    if (!{
      is_mm_dd_format || is_one_of_keys
    }) {
      cli_abort(paste(
        "If {.code highest_imputation = \"M\"} is specified, {.arg date_imputation} must be",
        "one of {.val first}, {.val mid}, {.val last}",
        "or a format with month and day specified as {.val mm-dd}: e.g. {.val 06-15}"
      ))
    }
  }

  if (highest_imputation == "Y") {
    assert_character_scalar(date_imputation, values = c("first", "last"))
  }
  date_imputation
}

#' Assert `time_imputation`
#'
#' Applies assertions on the `time_imputation` argument
#'
#' @param highest_imputation Highest imputation level
#' @param time_imputation The value to impute time when missing
#'
#' @returns asserted `time_imputation`
#'
#' @examples
#' # Assert valid 'first' time imputation
#' time_imp_first <- admiral:::assert_time_imputation("Y", "first")
#' print(time_imp_first)
#'
#' # Assert valid 'last' time imputation
#' time_imp_last <- admiral:::assert_time_imputation("Y", "last")
#' print(time_imp_last)
#'
#' # Assert valid custom time imputation "12:34:56"
#' time_imp_custom <- admiral:::assert_time_imputation("Y", "12:34:56")
#' print(time_imp_custom)
#'
#' @keywords internal
#'
#' @returns `invisible(NULL)`
assert_time_imputation <- function(highest_imputation, time_imputation) {
  time_imputation <-
    assert_character_scalar(
      time_imputation,
      case_sensitive = FALSE
    )

  time_imputation <- tolower(time_imputation)

  is_hh_mm_ss_format <- grepl(
    "^(0[0-9]|1[0-9]|2[0-3]|24):([0-5][0-9]):([0-5][0-9])$",
    time_imputation
  )
  is_one_of_keys <- time_imputation %in% c("first", "last")

  if (!{
    is_hh_mm_ss_format || is_one_of_keys
  } && highest_imputation != "n") {
    cli_abort(paste(
      "{.arg time_imputation} must be",
      "one of {.val first}, {.val last}",
      "or time specified as {.val hh:mm:ss}: e.g. {.val 12:00:00}"
    ))
  }

  return(time_imputation)
}

#' Assert Highest Imputation Validity
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function checks the validity and requirements for the `highest_imputation` parameter.
#' It ensures that necessary conditions are met when `highest_imputation` is set to "Y".
#'
#' @param highest_imputation A character scalar indicating the highest level of imputation.
#' @param highest_imputation_values A character vector of valid values for `highest_imputation`.
#' @param date_imputation Optional character scalar specifying the imputation method for dates.
#' @param max_dates Optional vector specifying maximum dates for imputation.
#' @param min_dates Optional vector specifying minimum dates for imputation.
#'
#' @details
#' - If `highest_imputation` is "Y", either `min_dates` or `max_dates` must be specified.
#' - If `highest_imputation` is "Y" and `date_imputation` is "first",
#' `min_dates` must be specified.
#' - If `highest_imputation` is "Y" and `date_imputation` is "last",
#' `max_dates` must be specified.
#'
#' @returns Returns `NULL` invisibly if assertions pass.
#'
#' @keywords internal
assert_highest_imputation <- function(highest_imputation, highest_imputation_values,
                                      date_imputation = NULL,
                                      max_dates, min_dates) {
  assert_character_scalar(
    highest_imputation,
    values = highest_imputation_values,
    case_sensitive = TRUE # not sure
  )

  if (highest_imputation != "Y") {
    return(invisible(NULL))
  }

  assert_character_scalar(date_imputation, values = c("first", "last"))
  no_mindates <- is.null(min_dates) || length(min_dates) == 0
  no_maxdates <- is.null(max_dates) || length(max_dates) == 0

  if (no_maxdates && no_mindates) {
    cli_abort(paste(
      "If {.code highest_imputation = \"Y\"} is specified, {.arg min_dates} or",
      "{.arg max_dates} must be specified respectively."
    ))
  }


  if (no_mindates && date_imputation == "first") {
    cli_abort(paste(
      "If {.code highest_imputation = \"Y\"} and {.code date_imputation = \"first\"}",
      "is specified, {.arg min_dates} must be specified."
    ))
  }

  if (no_maxdates && date_imputation == "last") {
    cli_abort(paste(
      "If {.code highest_imputation = \"Y\"} and {.code date_imputation = \"last\"}",
      "is specified, {.arg max_dates} must be specified."
    ))
  }


  return(invisible(NULL))
}

#' Get range of possible date / datetime
#'
#' @description
#' Internal helper function to convert a character vector of (possibly partial) dates (`dtc`)
#' into complete dates based on a specified imputation rule (`date_imputation`).
#'
#' @param dtc A character vector of dates in ISO 8601 format
#' (e.g., `"2022-12-15"`, `"2022-12"`, `"2022"`).
#' Partial dates are allowed.
#'
#' @param date_imputation A string specifying the imputation strategy to apply to
#' incomplete dates.
#' Accepts `"first"` (default) or `"last"` to impute to the first or last possible date,
#' respectively.
#'
#' @returns A character vector of fully imputed dates in `"YYYY-MM-DD"` format.
#'
#' @examples
#' # Impute date range for dates with 'first' date imputation
#' dtc_dates <- c("2020-02-29", "2021-03")
#' imputed_dates_first <- get_dt_dtm_range(dtc_dates, date_imputation = "first")
#' print(imputed_dates_first)
#'
#' # Impute date range for dates with 'last' date imputation
#' imputed_dates_last <- get_dt_dtm_range(dtc_dates, date_imputation = "last")
#' print(imputed_dates_last)
#'
#' # Impute datetime range with 'first' time imputation
#' dtc_datetimes <- c("2020-02-29T12:00", "2021-03T14:30")
#' imputed_datetimes_first <- get_dt_dtm_range(dtc_datetimes,
#'   date_imputation = "first",
#'   time_imputation = "first"
#' )
#' print(imputed_datetimes_first)
#'
#' # Impute datetime range with 'last' time imputation
#' imputed_datetimes_last <- get_dt_dtm_range(dtc_datetimes,
#'   date_imputation = "first",
#'   time_imputation = "last"
#' )
#' print(imputed_datetimes_last)
#'
#' # Edge case: Return empty character vector for empty input
#' imputed_empty <- get_dt_dtm_range(character(0), date_imputation = "first")
#' print(imputed_empty)
#'
#' @details
#' This function is a simplified version of the `impute_dtc_dt()` function that
#' does not throw an error when having no min_dates nor max_dates specified when
#' `date_imputation = "Y"`.
#' The function performs the following steps:
#' - Validates and parses the input date strings.
#' - Determines which components (year, month, day) are missing.
#' - Applies the selected imputation rule (`"first"` or `"last"`) to fill in
#'  missing components.
#' - Returns a fully qualified date or `NA` if imputation cannot be performed.
#'
#' @keywords internal
get_dt_dtm_range <- function(dtc,
                             date_imputation = "first",
                             time_imputation = NULL) {
  assert_character_vector(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  is_datetime <- !is.null(time_imputation)

  highest_imputation_level <- get_highest_imputation_level(is_datetime, "Y")

  date_imputation <- assert_character_scalar(date_imputation, case_sensitive = FALSE)
  if (is_datetime) {
    time_imputation <- assert_character_scalar(time_imputation, case_sensitive = FALSE)
    assert_time_imputation(
      highest_imputation = highest_imputation_level,
      time_imputation = time_imputation
    )
  }

  if (length(dtc) == 0) {
    return(character(0))
  }

  # Parse partials
  partial <- parse_partial_date_time(dtc, is_datetime)
  components <- names(partial)

  partial <- propagate_na_values(partial, is_datetime)

  target <- get_imputation_targets(partial, date_imputation, time_imputation, is_datetime)

  imputed <- impute_values(partial, target, components)
  imputed_dtc <- format_imputed_dtc(imputed, is_datetime)

  if (date_imputation == "last") {
    imputed_dtc <- adjust_last_day_imputation(imputed_dtc, partial, is_datetime)
  }

  return(imputed_dtc)
}

#' Get Highest Imputation Level
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' Determines the highest imputation level based on whether it's a date or datetime.
#'
#' @param is_datetime A logical indicating whether the imputation is for a datetime.
#' @param highest_imputation A character indicating the highest imputation level.
#'
#' @returns An integer representing the highest imputation level.
#'
#' @examples
#' # Get highest imputation level for date
#' highest_level_date <- get_highest_imputation_level(
#'   is_datetime = FALSE,
#'   highest_imputation = "Y"
#' )
#' print(highest_level_date)
#'
#' # Get highest imputation level for datetime
#' highest_level_datetime <- get_highest_imputation_level(
#'   is_datetime = TRUE,
#'   highest_imputation = "Y"
#' )
#' print(highest_level_datetime)
#'
#' # Get highest imputation level for date with month level
#' highest_level_month_date <- get_highest_imputation_level(
#'   is_datetime = FALSE,
#'   highest_imputation = "M"
#' )
#' print(highest_level_month_date)
#'
#' # Get highest imputation level for datetime with hour level
#' highest_level_hour_datetime <- get_highest_imputation_level(
#'   is_datetime = TRUE,
#'   highest_imputation = "h"
#' )
#' print(highest_level_hour_datetime)
#'
#' @keywords internal
get_highest_imputation_level <- function(is_datetime, highest_imputation) {
  if (is_datetime) dtm_level(highest_imputation) else dt_level(highest_imputation)
}

#' Get Imputation Targets
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' Determines the imputation targets for date and time components.
#'
#' @param partial A list of partial date/time components.
#' @param date_imputation A character specifying the date imputation method.
#' @param time_imputation A character specifying the time imputation method.
#' @param is_datetime A logical indicating whether it's a datetime imputation.
#'
#' @returns A list of imputation targets for date and (if applicable) time components.
#'
#' @examples
#' # Get imputation targets for a date with 'first' date imputation
#' partial_date <- list(year = "2020", month = "03", day = NA_character_)
#' target_first_date <- get_imputation_targets(partial_date,
#'   date_imputation = "first",
#'   time_imputation = NULL,
#'   is_datetime = FALSE
#' )
#' print(target_first_date)
#'
#' # Get imputation targets for a datetime with 'first' date and time imputation
#' partial_datetime <- list(
#'   year = "2020",
#'   month = "03",
#'   day = NA_character_,
#'   hour = "12",
#'   minute = NA_character_,
#'   second = NA_character_
#' )
#' target_first_datetime <- get_imputation_targets(partial_datetime,
#'   date_imputation = "first",
#'   time_imputation = "first",
#'   is_datetime = TRUE
#' )
#' print(target_first_datetime)
#'
#' # Get imputation targets for a datetime with 'last' date and time imputation
#' target_last_datetime <- get_imputation_targets(partial_datetime,
#'   date_imputation = "last",
#'   time_imputation = "last",
#'   is_datetime = TRUE
#' )
#' print(target_last_datetime)
#'
#' # Get imputation targets for a date with custom date imputation '06-15'
#' target_custom_date <- get_imputation_targets(partial_date,
#'   date_imputation = "06-15",
#'   time_imputation = NULL,
#'   is_datetime = FALSE
#' )
#' print(target_custom_date)
#'
#' # Get imputation targets for a datetime with custom time imputation '12:34:56'
#' target_custom_time <- get_imputation_targets(partial_datetime,
#'   date_imputation = "first",
#'   time_imputation = "12:34:56",
#'   is_datetime = TRUE
#' )
#' print(target_custom_time)
#'
#' @keywords internal
get_imputation_targets <- function(partial, date_imputation, time_imputation, is_datetime) {
  target_date <- get_imputation_target_date(
    date_imputation = date_imputation,
    month = partial[["month"]]
  )

  if (is_datetime) {
    target_time <- get_imputation_target_time(time_imputation = time_imputation)
    return(c(target_date, target_time))
  }

  return(target_date)
}

#' Adjust Last Day Imputation
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' Adjusts the imputed date/datetime to the last day of the month when necessary.
#'
#' @param imputed_dtc A character vector of imputed date/datetime strings.
#' @param partial A list of partial date/time components.
#' @param is_datetime A logical indicating whether it's a datetime adjustment.
#'
#' @returns A character vector of adjusted date/datetime strings.
#'
#' @examples
#' # Adjust last day imputation for a date with an incomplete day
#' imputed_date <- "2021-03-01"
#' partial_date <- list(year = "2021", month = "03", day = NA_character_)
#' adjusted_date <- admiral:::adjust_last_day_imputation(imputed_date,
#'   partial_date,
#'   is_datetime = FALSE
#' )
#' print(adjusted_date)
#'
#' # Adjust last day imputation for a datetime with missing day
#' imputed_datetime <- "2021-03-01T00:00:00"
#' partial_datetime <- list(
#'   year = "2021", month = "03", day = NA_character_,
#'   hour = "00", minute = "00", second = "00"
#' )
#' adjusted_datetime <- admiral:::adjust_last_day_imputation(imputed_datetime,
#'   partial_datetime,
#'   is_datetime = TRUE
#' )
#' print(adjusted_datetime)
#'
#' # Adjust last day imputation for a date with known day
#' partial_date_known_day <- list(year = "2021", month = "03", day = "15")
#' adjusted_date_known_day <- admiral:::adjust_last_day_imputation(imputed_date,
#'   partial_date_known_day,
#'   is_datetime = FALSE
#' )
#' print(adjusted_date_known_day)
#'
#' # Adjust last day imputation for a datetime with known day
#' partial_datetime_known_day <- list(
#'   year = "2021", month = "03", day = "15",
#'   hour = "00", minute = "00", second = "00"
#' )
#' adjusted_datetime_known_day <- admiral:::adjust_last_day_imputation(imputed_datetime,
#'   partial_datetime_known_day,
#'   is_datetime = TRUE
#' )
#' print(adjusted_datetime_known_day)
#'
#' @details
#' This function is used when the date imputation method is set to "last" and
#' the day component is missing. It adjusts the date to the last day of the month.
#'
#' @keywords internal
#'
#' @importFrom lubridate ymd_hms ymd rollback
adjust_last_day_imputation <- function(imputed_dtc, partial, is_datetime) {
  if_else(
    is.na(partial[["day"]]),
    if (is_datetime) {
      strftime(
        rollback(ymd_hms(imputed_dtc) + months(1)),
        format = "%Y-%m-%dT%H:%M:%S",
        tz = "UTC"
      )
    } else {
      strftime(
        rollback(ymd(imputed_dtc) + months(1)),
        format = "%Y-%m-%d"
      )
    },
    imputed_dtc
  )
}

#' Impute Missing Values
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' Imputes missing values in partial date/time components using target values.
#'
#' @param partial A list of partial date/time components.
#' @param target A list of target values for imputation.
#' @param components A character vector of component names.
#'
#' @returns A list of imputed date/time components.
#'
#' @examples
#' # Impute missing values for date components
#' partial_date <- list(year = "2020", month = NA_character_, day = NA_character_)
#' target_date <- list(year = "2020", month = "01", day = "01")
#' components_date <- c("year", "month", "day")
#' imputed_date <- impute_values(partial_date, target_date, components_date)
#' print(imputed_date)
#'
#' # Impute missing values for datetime components
#' partial_datetime <- list(
#'   year = "2020", month = NA_character_, day = NA_character_,
#'   hour = "12", minute = NA_character_, second = NA_character_
#' )
#' target_datetime <- list(
#'   year = "2020", month = "01", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' components_datetime <- c("year", "month", "day", "hour", "minute", "second")
#' imputed_datetime <- impute_values(partial_datetime, target_datetime, components_datetime)
#' print(imputed_datetime)
#'
#' # Impute missing values when some components are already present
#' partial_mixed <- list(year = "2020", month = "06", day = NA_character_)
#' target_mixed <- list(year = "2020", month = "01", day = "01")
#' components_mixed <- c("year", "month", "day")
#' imputed_mixed <- impute_values(partial_mixed, target_mixed, components_mixed)
#' print(imputed_mixed)
#'
#' @keywords internal
impute_values <- function(partial, target, components) {
  imputed <- vector("list", length(components))
  names(imputed) <- components
  for (c in components) {
    imputed[[c]] <- if_else(is.na(partial[[c]]), target[[c]], partial[[c]])
  }
  imputed
}

#' Format Imputed Date/Datetime
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' Formats imputed date/datetime components into a string representation.
#'
#' @param imputed A list of imputed date/time components.
#' @param is_datetime A logical indicating whether it's a datetime format.
#'
#' @returns A character vector of formatted date/datetime strings.
#'
#' @examples
#' # Format imputed datetime components
#' imputed_datetime <- list(
#'   year = "2020", month = "01", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' formatted_datetime <- format_imputed_dtc(imputed_datetime, is_datetime = TRUE)
#' print(formatted_datetime)
#'
#' # Format imputed date components
#' imputed_date <- list(year = "2020", month = "01", day = "01")
#' formatted_date <- format_imputed_dtc(imputed_date, is_datetime = FALSE)
#' print(formatted_date)
#'
#' # Handle imputed datetime with missing parts (contains 'x')
#' # Expected: NA because 'x' is an undefined component
#' imputed_partial_datetime <- list(
#'   year = "2020", month = "xx", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' formatted_partial_datetime <- format_imputed_dtc(imputed_partial_datetime, is_datetime = TRUE)
#' print(formatted_partial_datetime)
#'
#' # Handle imputed date with missing parts (contains 'x')
#' # Expected: NA because 'x' is an undefined component
#' imputed_partial_date <- list(year = "2020", month = "xx", day = "01")
#' formatted_partial_date <- format_imputed_dtc(imputed_partial_date, is_datetime = FALSE)
#' print(formatted_partial_date)
#'
#' @details
#' The function formats the imputed components into "YYYY-MM-DD" for dates
#' and "YYYY-MM-DDTHH:MM:SS" for datetimes. It replaces any string containing
#' 'x' with NA.
#'
#' @keywords internal
format_imputed_dtc <- function(imputed, is_datetime) {
  if (is_datetime) {
    dtc <- paste0(
      paste(imputed[["year"]], imputed[["month"]], imputed[["day"]], sep = "-"),
      "T",
      paste(imputed[["hour"]], imputed[["minute"]], imputed[["second"]], sep = ":")
    )
  } else {
    dtc <- paste(imputed[["year"]], imputed[["month"]], imputed[["day"]], sep = "-")
  }
  if_else(str_detect(dtc, "x"), NA_character_, dtc)
}

#' Propagate NA Values for datetime values
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`.
#' Propagates NA values through date/time components.
#'
#' @param partial A list of partial date/time components.
#'
#' @returns A list of date/time components with propagated NA values.
#'
#' @examples
#' # Propagate NA values through datetime components
#' partial_datetime <- list(
#'   year = "2020", month = NA_character_, day = "01",
#'   hour = "12", minute = NA_character_, second = "34"
#' )
#' propagated_datetime <- propagate_na_values(partial_datetime)
#' print(propagated_datetime)
#'
#' # Propagate NA values for datetime with missing higher order components
#' partial_missing <- list(
#'   year = NA_character_, month = "01", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' propagated_missing <- propagate_na_values(partial_missing)
#' print(propagated_missing)
#'
#' @details
#' This function ensures that if a higher-order component (e.g., month) is NA,
#' all lower-order components (e.g., day, hour, etc.) are also set to NA.
#'
#' @keywords internal
propagate_na_values <- function(partial, is_datetime = TRUE) {
  comp_length <- ifelse(is_datetime, 6, 3)
  for (i in 2:comp_length) {
    partial[[i]] <- if_else(is.na(partial[[i - 1]]), NA_character_, partial[[i]])
  }
  partial
}

#' Parse Partial Date or Datetime
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' This function parses a vector of date or datetime strings into their component
#' parts.
#'
#' @param dtc A character vector of date or datetime strings to be parsed.
#' @param is_datetime A logical value indicating whether the input strings include
#' time information.
#'
#' @returns A list of character vectors, each representing a component of the date
#' or datetime.
#' - For dates, the components are "year", "month", and "day".
#' - For datetimes, the components also include "hour", "minute", and "second".
#'
#' @examples
#' # Parse partial datetime components
#' dtc_datetime <- "2020-03-15T12:34"
#' parsed_datetime <- parse_partial_date_time(dtc_datetime, is_datetime = TRUE)
#' print(parsed_datetime)
#'
#' # Parse partial date components
#' dtc_date <- "2020-03"
#' parsed_date <- parse_partial_date_time(dtc_date, is_datetime = FALSE)
#' print(parsed_date)
#'
#' # Parse partial datetime with missing components
#' dtc_partial_datetime <- "2020-03T12"
#' parsed_partial_datetime <- parse_partial_date_time(dtc_partial_datetime, is_datetime = TRUE)
#' print(parsed_partial_datetime)
#'
#' # Parse partial date with incomplete year
#' dtc_year <- "2020"
#' parsed_year <- parse_partial_date_time(dtc_year, is_datetime = FALSE)
#' print(parsed_year)
#'
#' @details
#' The function uses different parsing methods depending on whether the input is
#'  a date or a datetime:
#' - For dates, it calls `get_partialdate()`.
#' - For datetimes, it calls `get_partialdatetime()`.
#'
#' @keywords internal
#'
parse_partial_date_time <- function(dtc, is_datetime) {
  if (is_datetime) {
    get_partialdatetime(dtc)
  } else {
    get_partialdatetime(dtc)[c("year", "month", "day")]
  }
}
