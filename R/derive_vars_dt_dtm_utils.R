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
#'
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

  if (is_datetime) {
    partial <- propagate_na_values(partial)
  }

  target <- get_imputation_targets(partial, date_imputation, time_imputation, is_datetime)

  imputed <- impute_values(partial, target, components)
  imputed_dtc <- format_imputed_dtc(imputed, is_datetime)

  if (date_imputation == "last") {
    imputed_dtc <- adjust_last_day_imputation(imputed_dtc, partial, is_datetime)
  }


  return(imputed_dtc)
}

#' Extract Partial Date Components
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This is a helping function for `get_dt_dtm_range()`
#' This function extracts year, month, and day components from a vector of
#' date strings.
#'
#'
#' @param dtc A character vector of date strings to be parsed.
#'
#' @returns A list with three elements: "year", "month", and "day", each containing
#'   a character vector of the respective date components.
#'
#' @details
#' The function uses regular expressions to extract date components. It handles
#' partial dates and propagates NA values for missing higher-order components.
#'
#' @keywords internal
get_partialdate <- function(dtc) {
  two <- "(\\d{2}|-?)"
  partialdate <- str_match(dtc, paste0("(\\d{4}|-?)-?", two, "-?", two))
  components <- c("year", "month", "day")
  names_list <- setNames(vector("list", length(components)), components)

  for (i in seq_along(components)) {
    x <- partialdate[, i + 1]
    x[x %in% c("-", "")] <- NA_character_
    names_list[[i]] <- x
  }

  for (i in 2:3) {
    names_list[[i]] <- if_else(is.na(names_list[[i - 1]]), NA_character_, names_list[[i]])
  }

  return(names_list)
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
#' Propagate NA Values
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
#' @details
#' This function ensures that if a higher-order component (e.g., month) is NA,
#' all lower-order components (e.g., day, hour, etc.) are also set to NA.
#'
#' @keywords internal
propagate_na_values <- function(partial) {
  for (i in 2:6) {
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
    get_partialdate(dtc)
  }
}
