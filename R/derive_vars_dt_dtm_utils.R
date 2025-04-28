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
#' @return A datetime object
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
#' @return `invisible(NULL)`
assert_dt_dtm_inputs <- function(new_vars_prefix, max_dates, min_dates, # nolint: cyclocomp_linter
                                 flag_imputation, flag_imputation_values,
                                 highest_imputation, highest_imputation_values,
                                 date_imputation = NULL) {
  assert_character_scalar(new_vars_prefix)
  assert_vars(max_dates, optional = TRUE)
  assert_vars(min_dates, optional = TRUE)

  assert_character_scalar(
    highest_imputation,
    values = highest_imputation_values,
    case_sensitive = FALSE # not sure
  )

  assert_character_scalar(
    flag_imputation,
    values = flag_imputation_values,
    case_sensitive = FALSE
  )

  if ((highest_imputation == "Y" && is.null(min_dates) && is.null(max_dates)) ||
    (highest_imputation == "Y" && length(min_dates) == 0 && length(max_dates) == 0)) {
    cli_abort(paste(
      "If {.code highest_imputation = \"Y\"} is specified, {.arg min_dates} or",
      "{.arg max_dates} must be specified respectively."
    ))
  }

  if (highest_imputation == "Y") {
    assert_character_scalar(date_imputation, values = c("first", "last"))
  }

  if (highest_imputation == "Y" && is.null(min_dates) && date_imputation == "first") {
    cli_warn(paste(
      "If {.code highest_imputation = \"Y\"} and {.code date_imputation = \"first\"}",
      "is specified, {.arg min_dates} should be specified."
    ))
  }

  if (highest_imputation == "Y" && is.null(max_dates) && date_imputation == "last") {
    cli_warn(paste(
      "If {.code highest_imputation = \"Y\"} and {.code date_imputation = \"last\"}",
      "is specified, {.arg max_dates} should be specified."
    ))
  }

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
#' @keywords internal
#'
#' @return `invisible(NULL)`
assert_date_imputation <- function(highest_imputation, date_imputation) {
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
  return(invisible(NULL))
}

#' Assert `time_imputation`
#'
#' Applies assertions on the `time_imputation` argument
#'
#' @param highest_imputation Highest imputation level
#' @param time_imputation The value to impute time when missing
#'
#' @keywords internal
#'
#' @return `invisible(NULL)`
assert_time_imputation <- function(highest_imputation, time_imputation) {
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
      'one of `"first"`, `"last`"',
      'or time specified as {.val hh:mm:ss}: e.g. {.val 12:00:00}'
    ))
  }

  return(invisible(NULL))
}
