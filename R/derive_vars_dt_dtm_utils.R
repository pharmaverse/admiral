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
#' level_day <- admiral:::dtm_level("D")
#' print(level_day)
#'
#' # Create a dtm_level object with level "h" (hour)
#' level_hour <- admiral:::dtm_level("h")
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
#'  - For `date_imputation = "mid"` `"xxxx"`, `"06"`, `"30"` if `month` is `NA`.
#'   otherwise `"15"` returned.
#'  - For `date_imputation = "last"` `"9999"`, `"12"`, `"28"` are returned.
#'  - For `date_imputation = "<mm>-<dd>"` `"xxxx"`, `"<mm>"`, `"<dd>"` are returned.
#'
#'  `"xxxx"` indicates that the component is undefined. If an undefined
#'  component occurs in the imputed `--DTC` value, the imputed `--DTC` value is set to
#'  `NA_character_` in the imputation functions.
#'
#' @examples
#' # Get imputation target for "first"
#' target_first <- admiral:::get_imputation_target_date("first", month = NA)
#' print(target_first)
#'
#' # Get imputation target for "mid" with specified month
#' target_mid <- admiral:::get_imputation_target_date("mid", month = "03")
#' print(target_mid)
#'
#' # Get imputation target for "mid" with NA month
#' target_mid_na <- admiral:::get_imputation_target_date("mid", month = NA)
#' print(target_mid_na)
#'
#' # Get imputation target for "last"
#' target_last <- admiral:::get_imputation_target_date("last", month = NA)
#' print(target_last)
#'
#' # Get imputation target for custom date imputation "06-15"
#' target_custom <- admiral:::get_imputation_target_date("06-15", month = NA)
#' print(target_custom)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_imputation_target_date <- function(date_imputation,
                                       month) {
  if (date_imputation == "first") {
    list(
      year = "0000",
      month = "01",
      day = "01"
    )
  } else if (date_imputation == "mid") {
    list(
      year = "xxxx",
      month = "06",
      day = if_else(is.na(month), "30", "15")
    )
  } else if (date_imputation == "last") {
    list(
      year = "9999",
      month = "12",
      day = "28"
    )
  } else {
    list(
      year = "xxxx",
      month = str_sub(date_imputation, 1, 2),
      day = str_sub(date_imputation, 4, 5)
    )
  }
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
#' target_first_time <- admiral:::get_imputation_target_time("first")
#' print(target_first_time)
#'
#' # Get imputation target for "last" time
#' target_last_time <- admiral:::get_imputation_target_time("last")
#' print(target_last_time)
#'
#' # Get imputation target for custom time imputation "12-34-56"
#' target_custom_time <- admiral:::get_imputation_target_time("12-34-56")
#' print(target_custom_time)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso  [impute_dtc_dtm()]
get_imputation_target_time <- function(time_imputation) {
  if (time_imputation == "first") {
    list(hour = "00", minute = "00", second = "00")
  } else if (time_imputation == "last") {
    list(hour = "23", minute = "59", second = "59")
  } else {
    list(
      hour = str_sub(time_imputation, 1, 2),
      minute = str_sub(time_imputation, 4, 5),
      second = str_sub(time_imputation, 7, -1)
    )
  }
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
    dt
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

#' Parse `--DTC` variable and Determine Components
#'
#' @param dtc The `'--DTC'` date to parse
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
#'   valid value for any component.
#'
#' @param create_datetime logical scalar. If `TRUE` returns Datetime components.
#'   If `FALSE` returns Date components.
#'
#'
#' @returns A list of character vectors. The elements of the list are named
#'   "year", "month", "day", "hour", "minute", and "second". Missing components
#'   are set to `NA_character_`.
#'
#' @details The function can be replaced by the parttime parser once it is
#'   available.
#'
#' @examples
#' # Datetime
#' # Get partial datetime components for a complete datetime string
#' dtc_complete <- admiral:::get_partialdatetime("2020-03-15T12:34:56", TRUE)
#' print(dtc_complete)
#'
#' # Get partial datetime components for a partial datetime string
#' dtc_partial <- admiral:::get_partialdatetime("2020-03-15T12:34", TRUE)
#' print(dtc_partial)
#'
#' # Get partial datetime components for a date-only string
#' dtc_date_only <- admiral:::get_partialdatetime("2020-03-15", TRUE)
#' print(dtc_date_only)
#'
#' # Get partial datetime components for an incomplete year string
#' dtc_year_partial <- admiral:::get_partialdatetime("2020", TRUE)
#' print(dtc_year_partial)
#'
#' # Date
#' # Get partial date components for a complete datetime string
#' dtc_complete <- admiral:::get_partialdatetime("2020-03-15T12:34:56", FALSE)
#' print(dtc_complete)
#'
#' # Get partial date components for a partial datetime string
#' dtc_partial <- admiral:::get_partialdatetime("2020-03-15T12:34", FALSE)
#' print(dtc_partial)
#'
#' # Get partial date components for a year and month only string
#' dtc_month_only <- admiral:::get_partialdatetime("2020-03", FALSE)
#' print(dtc_month_only)
#'
#' # Get partial date components for an incomplete year string
#' dtc_year_partial <- admiral:::get_partialdatetime("2020", FALSE)
#' print(dtc_year_partial)
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_partialdatetime <- function(dtc, create_datetime) {
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

  if (!create_datetime) {
    partial <- partial[c("year", "month", "day")]
  }

  partial
}

#' Assert `date_imputation`
#'
#' @description
#' Applies assertions on the `date_imputation` argument to reduce
#' cyclomatic complexity
#'
#' @details
#' Asserts that date_imputation is a scalar.
#' Asserts that the values in `date_imputation` are permitted.
#' The permitted values in `date_imputation` vary by `highest_imputation`
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#' @param highest_imputation Highest imputation level
#'
#' @returns asserted `date_imputation`
#'
#' @keywords internal
assert_date_imputation <- function(date_imputation, highest_imputation) {
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
    is_mm_dd_format <- str_detect(date_imputation, "^(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$")
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
#' @param time_imputation The value to impute time when missing
#' @param highest_imputation Highest imputation level
#'
#' @returns asserted `time_imputation`
#'
#' @examples
#' # Assert valid 'first' time imputation
#' time_imp_first <- admiral:::assert_time_imputation("first", "Y")
#' print(time_imp_first)
#'
#' # Assert valid 'last' time imputation
#' time_imp_last <- admiral:::assert_time_imputation("last", "Y")
#' print(time_imp_last)
#'
#' # Assert valid custom time imputation "12:34:56"
#' time_imp_custom <- admiral:::assert_time_imputation("12:34:56", "Y")
#' print(time_imp_custom)
#'
#' @keywords internal
#'
assert_time_imputation <- function(time_imputation, highest_imputation) {
  time_imputation <-
    assert_character_scalar(
      time_imputation,
      case_sensitive = FALSE
    )

  time_imputation <- tolower(time_imputation)

  is_hh_mm_ss_format <- str_detect(
    time_imputation,
    "^((0[0-9]|1[0-9]|2[0-3]):([0-5][0-9]):([0-5][0-9])|24:00:00)$"
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

  time_imputation
}

#' Assert Highest Imputation Validity
#'
#' @description
#'
#' This function checks the validity and requirements for the `highest_imputation` argument.
#' It ensures that necessary conditions for `date_imputation`, `min_dates`,
#' and `max_dates` are met when `highest_imputation` is set to `"Y"`.
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
    case_sensitive = TRUE
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
  invisible(NULL)
}

#' Get Range of Partial Date / Datetime
#'
#' @description
#' Internal helper function to convert a character vector of (possibly partial) dates (`dtc`)
#' into complete dates based on a specified imputation rule (`date_imputation`).
#'
#' @param dtc A character vector of dates in ISO 8601 format
#' (e.g., `"2022-12-15"`, `"2022-12"`, `"2022"`).
#' Partial dates are allowed.
#'
#' @param create_datetime return the range in datetime format.
#'
#' @returns A list containing two vectors of fully imputed dates
#' in `"YYYY-MM-DD"` or `"YYYY-MM-DDThh:mm:ss"` format - the lower and upper limit of the range.
#'
#' @examples
#' # Get Range from Partial Dates
#' dtc_dates <- c("2020-02-29", "2021-03")
#' imputed_dates_first <- admiral:::get_dt_dtm_range(dtc_dates, create_datetime = FALSE)
#' print(imputed_dates_first)
#'
#'
#' # Get Range from Partial Datetime
#' dtc_datetimes <- c("2020-02-29T12:00", "2021-03T14:30")
#' imputed_datetimes_first <- admiral:::get_dt_dtm_range(dtc_datetimes, create_datetime = TRUE)
#' print(imputed_datetimes_first)
#'
#' # Edge case: Return empty character vector for empty input
#' imputed_empty <- admiral:::get_dt_dtm_range(character(0), create_datetime = TRUE)
#' print(imputed_empty)
#'
#' @details
#' The functions replaces missing components in `dtc` with the earliest (lower bound)
#' and latest (upper bound) possible value. Missing year is replaced with `"0000"` for the
#' lower bound and `"9999"` for the upper bound.
#'
#' @keywords internal
get_dt_dtm_range <- function(dtc,
                             create_datetime) {
  assert_character_vector(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  highest_imputation_level <- get_highest_imputation_level(
    highest_imputation = "Y",
    create_datetime = create_datetime
  )

  if (length(dtc) == 0) {
    return(character(0))
  }

  # Parse partials
  partial <- get_partialdatetime(dtc, create_datetime = create_datetime)
  partial <- propagate_na_values(partial)

  lo_up <- c("first", "last")

  targets <- lapply(lo_up, function(x) {
    get_imputation_targets(partial,
      date_imputation = x,
      time_imputation = x
    )
  })

  imputed_targets <- lapply(targets, function(target) {
    impute_date_time(partial, target)
  })

  imputed_dtcs <- lapply(imputed_targets, function(imputed) {
    format_imputed_dtc(imputed)
  })

  imputed_dtcs[[2]] <- adjust_last_day_imputation(imputed_dtcs[[2]], partial)

  names(imputed_dtcs) <- c("lower", "upper")

  imputed_dtcs
}

#' Get Highest Imputation Level
#'
#' @description
#'
#' Returns the `dt_level()` or `dtm_level()` representation of the `highest_imputation`
#' character value. The level object allows comparisons of levels.
#'
#' @param highest_imputation A character indicating the highest imputation level.
#' @param create_datetime A logical indicating whether datetime factors levels are required.
#'
#' @returns A `dt_level()` or `dtm_level()` object representing the highest imputation level.
#'
#' @examples
#' # Get highest imputation level for date
#' highest_level_date <- admiral:::get_highest_imputation_level(
#'   highest_imputation = "Y",
#'   create_datetime = FALSE
#' )
#' print(highest_level_date)
#'
#' # Get highest imputation level for datetime
#' highest_level_datetime <- admiral:::get_highest_imputation_level(
#'   highest_imputation = "Y",
#'   create_datetime = TRUE
#' )
#' print(highest_level_datetime)
#'
#' # Get highest imputation level for date with month level
#' highest_level_month_date <- admiral:::get_highest_imputation_level(
#'   highest_imputation = "M",
#'   create_datetime = FALSE
#' )
#' print(highest_level_month_date)
#'
#' # Get highest imputation level for datetime with hour level
#' highest_level_hour_datetime <- admiral:::get_highest_imputation_level(
#'   highest_imputation = "h",
#'   create_datetime = TRUE
#' )
#' print(highest_level_hour_datetime)
#'
#' @keywords internal
get_highest_imputation_level <- function(highest_imputation, create_datetime) {
  if (create_datetime) dtm_level(highest_imputation) else dt_level(highest_imputation)
}

#' Get Imputation Targets
#'
#' @description
#'
#' Determines the imputation targets for date (see `get_imputation_target_date()` and time
#' (see `get_imputation_target_time()`) components.
#'
#' @param partial A list of partial date/time components.
#' @inheritParams get_imputation_target_date
#' @inheritParams get_imputation_target_time
#'
#' @returns A list of imputation targets for date and (if applicable) time components.
#'
#' @examples
#' # Get imputation targets for a date with 'first' date imputation
#' partial_date <- list(year = "2020", month = "03", day = NA_character_)
#' target_first_date <- admiral:::get_imputation_targets(partial_date,
#'   date_imputation = "first",
#'   time_imputation = NULL
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
#' target_first_datetime <- admiral:::get_imputation_targets(partial_datetime,
#'   date_imputation = "first",
#'   time_imputation = "first"
#' )
#' print(target_first_datetime)
#'
#' # Get imputation targets for a datetime with 'last' date and time imputation
#' target_last_datetime <- admiral:::get_imputation_targets(partial_datetime,
#'   date_imputation = "last",
#'   time_imputation = "last"
#' )
#' print(target_last_datetime)
#'
#' # Get imputation targets for a date with custom date imputation '06-15'
#' target_custom_date <- admiral:::get_imputation_targets(partial_date,
#'   date_imputation = "06-15",
#'   time_imputation = NULL
#' )
#' print(target_custom_date)
#'
#' # Get imputation targets for a datetime with custom time imputation '12:34:56'
#' target_custom_time <- admiral:::get_imputation_targets(partial_datetime,
#'   date_imputation = "first",
#'   time_imputation = "12:34:56"
#' )
#' print(target_custom_time)
#'
#' @keywords internal
get_imputation_targets <- function(partial, date_imputation = NULL, time_imputation = NULL) {
  target_date <- get_imputation_target_date(
    date_imputation = date_imputation,
    month = partial[["month"]]
  )

  if (is.null(time_imputation) && is_partial_datetime(partial)) {
    cli_abort("As {.arg partial} is datetime, {.arg time_imputation} is expected.")
  }
  if (is_partial_datetime(partial)) {
    target_time <- get_imputation_target_time(time_imputation = time_imputation)
    return(c(target_date, target_time))
  }

  target_date
}

#' Adjust Last Day Imputation
#'
#' @description
#'
#' This functions adjusts the day of the imputed date to the last day the month
#' if the day was imputed. It should be called if `date_imputation = "last"` was used
#' for the date imputation as `get_imputation_target_date()` imputes the last day
#' as `"28"`.
#'
#' @param imputed_dtc A character vector of imputed date/datetime strings.
#' @param partial A list of partial date/time components.
#'
#' @returns A character vector of adjusted date/datetime strings.
#'
#' @examples
#' # Adjust last day imputation for a date with an incomplete day
#' imputed_date <- "2021-03-28"
#' partial_date <- list(year = "2021", month = "03", day = NA_character_)
#' admiral:::adjust_last_day_imputation(imputed_date, partial_date)
#'
#' # Adjust last day imputation for a datetime with missing day
#' imputed_datetime <- "2021-03-28T00:00:00"
#' partial_datetime <- list(
#'   year = "2021", month = "03", day = NA_character_,
#'   hour = "00", minute = "00", second = "00"
#' )
#' admiral:::adjust_last_day_imputation(imputed_datetime, partial_datetime)
#'
#' # Adjust last day imputation for a date with known day
#' partial_date_known_day <- list(year = "2021", month = "03", day = "15")
#' adjusted_date_known_day <- admiral:::adjust_last_day_imputation(
#'   imputed_date,
#'   partial_date_known_day
#' )
#' print(adjusted_date_known_day)
#'
#' # Adjust last day imputation for a datetime with known day
#' partial_datetime_known_day <- list(
#'   year = "2021", month = "03", day = "15",
#'   hour = "00", minute = "00", second = "00"
#' )
#' adjusted_datetime_known_day <- admiral:::adjust_last_day_imputation(
#'   imputed_datetime,
#'   partial_datetime_known_day
#' )
#' print(adjusted_datetime_known_day)
#'
#' @details
#' If the day component in `partial` is missing,
#' the day (in `imputed_dtc`) is adjusted to the last day of the month.
#'
#' @keywords internal
#'
adjust_last_day_imputation <- function(imputed_dtc, partial) {
  if_else(
    is.na(partial[["day"]]),
    if (is_partial_datetime(partial)) {
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
#'
#' Imputes missing values in partial date/time components using target values.
#'
#' @param partial A list of partial date/time components.
#' @param target A list of target values for imputation.
#'
#' @returns A list of imputed date/time components.
#'
#' @examples
#' # Impute missing values for date components
#' partial_date <- list(year = "2020", month = NA_character_, day = NA_character_)
#' target_date <- list(year = "2020", month = "01", day = "01")
#' imputed_date <- admiral:::impute_date_time(partial_date, target_date)
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
#' imputed_datetime <- admiral:::impute_date_time(
#'   partial_datetime, target_datetime
#' )
#' print(imputed_datetime)
#'
#' # Impute missing values when some components are already present
#' partial_mixed <- list(year = "2020", month = "06", day = NA_character_)
#' target_mixed <- list(year = "2020", month = "01", day = "01")
#' imputed_mixed <- admiral:::impute_date_time(partial_mixed, target_mixed)
#' print(imputed_mixed)
#'
#' @keywords internal
impute_date_time <- function(partial, target) {
  # assert partial and target are consistent
  if (!isTRUE(all.equal(names(partial), names(target)))) {
    cli_abort("Names of {.arg partial} and {.arg target} do not match.")
  }
  components <- names(partial)
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
#'
#' Formats imputed date/datetime components into a string representation.
#'
#' @param imputed A list of imputed date/time components.
#'
#' @returns A character vector of formatted date/datetime strings.
#'
#' @examples
#' # Format imputed datetime components
#' imputed_datetime <- list(
#'   year = "2020", month = "01", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' formatted_datetime <- admiral:::format_imputed_dtc(imputed_datetime)
#' print(formatted_datetime)
#'
#' # Format imputed date components
#' imputed_date <- list(year = "2020", month = "01", day = "01")
#' formatted_date <- admiral:::format_imputed_dtc(imputed_date)
#' print(formatted_date)
#'
#' # Handle imputed datetime with missing parts (contains 'x')
#' # Expected: NA because 'x' is an undefined component
#' imputed_partial_datetime <- list(
#'   year = "2020", month = "xx", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' formatted_partial_datetime <- admiral:::format_imputed_dtc(imputed_partial_datetime)
#' print(formatted_partial_datetime)
#'
#' # Handle imputed date with missing parts (contains 'x')
#' # Expected: NA because 'x' is an undefined component
#' imputed_partial_date <- list(year = "2020", month = "xx", day = "01")
#' formatted_partial_date <- admiral:::format_imputed_dtc(imputed_partial_date)
#' print(formatted_partial_date)
#'
#' @details
#' The function formats the imputed components into `"YYYY-MM-DD"` for dates
#' and `"YYYY-MM-DDThh:mm:ss"` for datetimes. It replaces any string containing
#' `"x"` with `NA`.
#'
#' @keywords internal
format_imputed_dtc <- function(imputed) {
  if (is_partial_datetime(imputed)) {
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
#'
#' Propagates `NA` values through date/time components.
#'
#' @param partial A list of partial date/time components.
#'
#' @returns A list of date/time components with propagated `NA` values.
#'
#' @examples
#' # Propagate NA values through datetime components
#' partial_datetime <- list(
#'   year = "2020", month = NA_character_, day = "01",
#'   hour = "12", minute = NA_character_, second = "34"
#' )
#' propagated_datetime <- admiral:::propagate_na_values(partial_datetime)
#' print(propagated_datetime)
#'
#' # Propagate NA values for datetime with missing higher order components
#' partial_missing <- list(
#'   year = NA_character_, month = "01", day = "01",
#'   hour = "12", minute = "00", second = "00"
#' )
#' propagated_missing <- admiral:::propagate_na_values(partial_missing)
#' print(propagated_missing)
#'
#' partial_missing_date <- list(
#'   year = "2023", month = NA_character_, day = "01"
#' )
#' propagated_missing_date <- admiral:::propagate_na_values(partial_missing_date)
#' print(propagated_missing_date)
#'
#' @details
#' This function ensures that if a higher-order component (e.g., month) is `NA`,
#' all lower-order components (e.g., day, hour, etc.) are also set to `NA`.
#'
#' @keywords internal
propagate_na_values <- function(partial) {
  for (i in 2:length(partial)) {
    partial[[i]] <- if_else(is.na(partial[[i - 1]]), NA_character_, partial[[i]])
  }
  partial
}

#' Check if a Partial Date/Time is a Datetime
#'
#' @description
#' This function determines whether a given partial date/time structure represents
#' a datetime or just a date.
#'
#' @param partial A named list containing date or datetime components.
#'
#' @return A logical value. TRUE if the partial represents a datetime,
#' FALSE if it represents a date only.
#'
#' @details
#' The function checks for the presence of all date components (year, month, day)
#' and all time components (hour, minute, second) in the input list. If all components
#'  are present, it's considered a datetime.
#' If only date components are present, it's considered a date.
#' Any other combination will result in an error.
#'
#' @examples
#' # Datetime example
#' partial_datetime <- list(
#'   year = "2023", month = "05", day = "15",
#'   hour = "14", minute = "30", second = "00"
#' )
#' admiral:::is_partial_datetime(partial_datetime) # Returns TRUE
#'
#' # Date example
#' partial_date <- list(year = "2023", month = "05", day = "15")
#' admiral:::is_partial_datetime(partial_date) # Returns FALSE
#'
#' # Invalid example
#' \dontrun{
#' partial_invalid <- list(year = "2023", month = "05", hour = "14")
#' admiral:::is_partial_datetime(partial_invalid) # Throws an error
#' }
#'
#' @keywords internal
is_partial_datetime <- function(partial) {
  date_components <- c("year", "month", "day")
  time_components <- c("hour", "minute", "second")

  if (all(c(date_components, time_components) %in% names(partial))) {
    TRUE
  } else if (all(date_components %in% names(partial))) {
    FALSE
  } else {
    cli_abort(paste0(
      "{.arg partial} must be a named list containing either all date components",
      " or all datetime components"
    ))
  }
}
