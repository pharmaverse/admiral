#' Impute Partial Date(-time) Portion of a `'--DTC'` Variable
#'
#' Imputation partial date/time portion of a `'--DTC'` variable. based on user
#' input.
#'
#' @param dtc The `'--DTC'` date to impute
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
#'   valid "missing" value for any component.
#'
#' @param highest_imputation Highest imputation level
#'
#'   The `highest_imputation` argument controls which components of the DTC
#'   value are imputed if they are missing. All components up to the specified
#'   level are imputed.
#'
#'   If a component at a higher level than the highest imputation level is
#'   missing, `NA_character_` is returned. For example, for `highest_imputation
#'   = "D"` `"2020"` results in `NA_character_` because the month is missing.
#'
#'   If `"n"` is specified, no imputation is performed, i.e., if any component is
#'   missing, `NA_character_` is returned.
#'
#'   If `"Y"` is specified, `date_imputation` should be `"first"` or `"last"`
#'   and `min_dates` or `max_dates` should be specified respectively. Otherwise,
#'   `NA_character_` is returned if the year component is missing.
#'
#'   *Default*: `"h"`
#'
#'   *Permitted Values*: `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"h"` (hour), `"m"` (minute), `"s"` (second), `"n"` (none, lowest
#'   level)
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for the
#'   15th of June (The year can not be specified; for imputing the year
#'   `"first"` or `"last"` together with `min_dates` or `max_dates` argument can
#'   be used (see examples).),
#'   - or as a keyword: `"first"`, `"mid"`, `"last"` to impute to the first/mid/last
#'   day/month.
#'
#'   The argument is ignored if `highest_imputation` is less then `"D"`.
#'
#'   *Default*: `"first"`.
#'
#' @param time_imputation The value to impute the time when a timepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as `"hh:mm:ss"`: e.g. `"00:00:00"`
#'   for the start of the day,
#'   - or as a keyword: `"first"`,`"last"` to impute to the start/end of a day.
#'
#'   The argument is ignored if `highest_imputation = "n"`.
#'
#'   *Default*: `"first"`.
#'
#' @param min_dates Minimum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not
#' before any of the specified dates, e.g., that the imputed adverse event start
#' date is not before the first treatment date. Only dates which are in the
#' range of possible dates of the `dtc` value are considered. The possible dates
#' are defined by the missing parts of the `dtc` date (see example below). This
#' ensures that the non-missing parts of the `dtc` date are not changed.
#' A date or date-time object is expected.
#' For example
#'
#' ```{r echo=TRUE, eval=FALSE}
#' impute_dtc_dtm(
#'   "2020-11",
#'   min_dates = list(
#'    ymd_hms("2020-12-06T12:12:12"),
#'    ymd_hms("2020-11-11T11:11:11")
#'   ),
#'   highest_imputation = "M"
#' )
#' ```
#'
#' returns `"2020-11-11T11:11:11"` because the possible dates for `"2020-11"`
#' range from `"2020-11-01T00:00:00"` to `"2020-11-30T23:59:59"`. Therefore
#' `"2020-12-06T12:12:12"` is ignored. Returning `"2020-12-06T12:12:12"` would
#' have changed the month although it is not missing (in the `dtc` date).
#'
#' For date variables (not datetime) in the list the time is imputed to
#' `"00:00:00"`. Specifying date variables makes sense only if the date is
#' imputed. If only time is imputed, date variables do not affect the result.
#'
#' @param max_dates Maximum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not after
#' any of the specified dates, e.g., that the imputed date is not after the data
#' cut off date. Only dates which are in the range of possible dates are
#' considered. A date or date-time object is expected.
#'
#' For date variables (not datetime) in the list the time is imputed to
#' `"23:59:59"`. Specifying date variables makes sense only if the date is
#' imputed. If only time is imputed, date variables do not affect the result.

#' @param preserve Preserve day if month is missing and day is present
#'
#' For example `"2019---07"` would return `"2019-06-07` if `preserve = TRUE`
#' (and `date_imputation = "mid"`).
#'
#' Permitted Values: `TRUE`, `FALSE`
#'
#' *Default*: `FALSE`
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return A character vector
#'
#' @author Samia Kabi, Stefan Bundfuss
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' dates <- c(
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18T15",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' # No date imputation (highest_imputation defaulted to "h")
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc_dtm(dtc = dates)
#'
#' # No date imputation (highest_imputation defaulted to "h")
#' # Missing time part imputed with 23:59:59 portion
#' impute_dtc_dtm(
#'   dtc = dates,
#'   time_imputation = "23:59:59"
#' )
#'
#' # Same as above
#' impute_dtc_dtm(
#'   dtc = dates,
#'   time_imputation = "last"
#' )
#'
#' # Impute to first day/month if date is partial
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc_dtm(
#'   dtc = dates,
#'   highest_imputation = "M"
#' )
#' # same as above
#' impute_dtc_dtm(
#'   dtc = dates,
#'   highest_imputation = "M",
#'   date_imputation = "01-01"
#' )
#'
#' # Impute to last day/month if date is partial
#' # Missing time part imputed with 23:59:59 portion
#' impute_dtc_dtm(
#'   dtc = dates,
#'   date_imputation = "last",
#'   time_imputation = "last"
#' )
#'
#' # Impute to mid day/month if date is partial
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc_dtm(
#'   dtc = dates,
#'   highest_imputation = "M",
#'   date_imputation = "mid"
#' )
#'
#' # Impute a date and ensure that the imputed date is not before a list of
#' # minimum dates
#' impute_dtc_dtm(
#'   "2020-12",
#'   min_dates = list(
#'     ymd_hms("2020-12-06T12:12:12"),
#'     ymd_hms("2020-11-11T11:11:11")
#'   ),
#'   highest_imputation = "M"
#' )
#'
#' # Impute completely missing dates (only possible if min_dates or max_dates is specified)
#' impute_dtc_dtm(
#'   c("2020-12", NA_character_),
#'   min_dates = list(
#'     ymd_hms("2020-12-06T12:12:12"),
#'     ymd_hms("2020-11-11T11:11:11")
#'   ),
#'   highest_imputation = "Y"
#' )
impute_dtc_dtm <- function(dtc,
                           highest_imputation = "h",
                           date_imputation = "first",
                           time_imputation = "first",
                           min_dates = NULL,
                           max_dates = NULL,
                           preserve = FALSE) {
  # Check arguments ----
  assert_character_vector(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)
  imputation_levels <- c(
    none = "n",
    second = "s",
    minute = "m",
    hour = "h",
    day = "D",
    month = "M",
    year = "Y"
  )
  assert_character_scalar(highest_imputation, values = imputation_levels)
  highest_imputation <- dtm_level(highest_imputation)
  date_imputation <-
    assert_character_scalar(
      date_imputation,
      case_sensitive = FALSE
    )
  time_imputation <-
    assert_character_scalar(
      time_imputation,
      case_sensitive = FALSE
    )
  assert_logical_scalar(preserve)

  if (length(dtc) == 0) {
    return(vector("character"))
  }

  # Parse character date ----
  partial <- get_partialdatetime(dtc)
  components <- names(partial)

  # Handle preserve argument ----
  if (!preserve) {
    for (i in 2:6) {
      partial[[i]] <- if_else(is.na(partial[[i - 1]]), NA_character_, partial[[i]])
    }
  }
  # Determine target components ----
  target_date <- get_imputation_target_date(
    date_imputation = date_imputation,
    month = partial[["month"]]
  )
  target_time <- get_imputation_target_time(
    time_imputation = time_imputation
  )
  target <- c(target_date, target_time)

  for (c in components) {
    if (highest_imputation < dtm_level(imputation_levels[[c]])) {
      target[[c]] <- "xx"
    }
  }

  # Impute ----
  imputed <- vector("list", 6)
  names(imputed) <- components
  for (c in components) {
    imputed[[c]] <- if_else(is.na(partial[[c]]), target[[c]], partial[[c]])
  }

  imputed_dtc <-
    paste0(
      paste(imputed[["year"]], imputed[["month"]], imputed[["day"]], sep = "-"),
      "T",
      paste(imputed[["hour"]], imputed[["minute"]], imputed[["second"]], sep = ":")
    )

  imputed_dtc <-
    if_else(
      str_detect(imputed_dtc, "x"),
      NA_character_,
      imputed_dtc
    )

  if (date_imputation == "last") {
    imputed_dtc <-
      if_else(
        is.na(partial[["day"]]),
        strftime(
          rollback(ymd_hms(imputed_dtc) + months(1)),
          format = "%Y-%m-%dT%H:%M:%S",
          tz = "UTC"
        ),
        imputed_dtc
      )
  }

  # Handle min_dates and max_dates argument ----
  restrict_imputed_dtc_dtm(
    dtc,
    imputed_dtc = imputed_dtc,
    min_dates = min_dates,
    max_dates = max_dates
  )
}

#' Create a `dtm_level` object
#'
#' @param level Datetime level
#'
#'   *Permitted Values*: `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"h"` (hour), `"m"` (minute), `"s"` (second, lowest level), `"n"`
#'   (none)
#'
#' @returns A `dtm_level` object
#'
#' @details A `dtm_level` object is an ordered factor, i.e., two objects can be
#'   compared.
#'
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
dtm_level <- function(level) {
  out <-
    factor(
      level,
      levels = c("n", "s", "m", "h", "D", "M", "Y"),
      ordered = TRUE
    )
  class(out) <- c("dtm_level", class(out))
  out
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
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_partialdatetime <- function(dtc) {
  two <- "(\\d{2}|-?)"
  partialdate <- stringr::str_match(dtc, paste0(
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
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
get_imputation_target_date <- function(date_imputation,
                                       month) {
  target <- vector("list", 3)
  names(target) <- c("year", "month", "day")
  if (date_imputation == "first") {
    target[["year"]] <- "0000"
    target[["month"]] <- "01"
    target[["day"]] <- "01"
  } else if (date_imputation == "mid") {
    target[["year"]] <- "xxxx"
    target[["month"]] <- "06"
    target[["day"]] <- if_else(is.na(month), "30", "15")
  } else if (date_imputation == "last") {
    target[["year"]] <- "9999"
    target[["month"]] <- "12"
    target[["day"]] <- "28"
  } else {
    target[["year"]] <- "xxxx"
    target[["month"]] <- str_sub(date_imputation, 1, 2)
    target[["day"]] <- str_sub(date_imputation, 4, 5)
  }
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
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
#'
#' @seealso  [impute_dtc_dtm()]
get_imputation_target_time <- function(time_imputation) {
  target <- vector("list", 3)
  names(target) <- c("hour", "minute", "second")
  if (time_imputation == "first") {
    target[["hour"]] <- "00"
    target[["minute"]] <- "00"
    target[["second"]] <- "00"
  } else if (time_imputation == "last") {
    target[["hour"]] <- "23"
    target[["minute"]] <- "59"
    target[["second"]] <- "59"
  } else {
    target[["hour"]] <- str_sub(time_imputation, 1, 2)
    target[["minute"]] <- str_sub(time_imputation, 4, 5)
    target[["second"]] <- str_sub(time_imputation, 7, -1)
  }
  target
}

#' Restrict Imputed DTC date to Minimum/Maximum Dates
#'
#' @param imputed_dtc The imputed DTC date
#'
#' @inheritParams impute_dtc_dtm
#'
#' @returns
#'   - The last of the minimum dates (`min_dates`) which are in the range of the
#'   partial DTC date (`dtc`)
#'   - The first of the maximum dates (`max_dates`) which are in the range of the
#'   partial DTC date (`dtc`)
#'   - `imputed_dtc` if the partial DTC date (`dtc`) is not in range of any of
#'   the minimum or maximum dates.
#'
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
restrict_imputed_dtc_dtm <- function(dtc,
                                     imputed_dtc,
                                     min_dates,
                                     max_dates) {
  if (!is.null(min_dates) | !is.null(max_dates)) {
    # determine range of possible dates
    min_dtc <-
      impute_dtc_dtm(
        dtc,
        highest_imputation = "Y",
        date_imputation = "first",
        time_imputation = "first"
      )
    max_dtc <-
      impute_dtc_dtm(
        dtc,
        highest_imputation = "Y",
        date_imputation = "last",
        time_imputation = "last"
      )
  }
  if (!is.null(min_dates)) {
    # for each minimum date within the range ensure that the imputed date is not
    # before it
    for (min_date in min_dates) {
      assert_date_vector(min_date)
      min_date_iso <- strftime(min_date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
      imputed_dtc <- if_else(
        min_dtc <= min_date_iso & min_date_iso <= max_dtc,
        pmax(imputed_dtc, min_date_iso),
        imputed_dtc,
        missing = imputed_dtc
      )
    }
  }
  if (!is.null(max_dates)) {
    # for each maximum date within the range ensure that the imputed date is not
    # after it
    for (max_date in max_dates) {
      assert_date_vector(max_date)
      max_date <- convert_date_to_dtm(
        max_date,
        time_imputation = "last"
      )
      max_date_iso <- strftime(max_date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
      imputed_dtc <- if_else(
        min_dtc <= max_date_iso & max_date_iso <= max_dtc,
        pmin(imputed_dtc, max_date_iso),
        imputed_dtc,
        missing = imputed_dtc
      )
    }
  }
  imputed_dtc
}

#' Impute Partial Date Portion of a `'--DTC'` Variable
#'
#' Imputation partial date portion of a `'--DTC'` variable based on user input.
#'
#' @param dtc The `'--DTC'` date to impute
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
#'   valid "missing" value for any component.
#'
#' @param highest_imputation Highest imputation level
#'
#'   The `highest_imputation` argument controls which components of the DTC
#'   value are imputed if they are missing. All components up to the specified
#'   level are imputed.
#'
#'   If a component at a higher level than the highest imputation level is
#'   missing, `NA_character_` is returned. For example, for `highest_imputation
#'   = "D"` `"2020"` results in `NA_character_` because the month is missing.
#'
#'   If `"n"` is specified no imputation is performed, i.e., if any component is
#'   missing, `NA_character_` is returned.
#'
#'   If `"Y"` is specified, `date_imputation` should be `"first"` or `"last"`
#'   and `min_dates` or `max_dates` should be specified respectively. Otherwise,
#'   `NA_character_` is returned if the year component is missing.
#'
#'   *Default*: `"n"`
#'
#'   *Permitted Values*: `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"n"` (none, lowest level)
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for the
#'   15th of June (The year can not be specified; for imputing the year
#'   `"first"` or `"last"` together with `min_dates` or `max_dates` argument can
#'   be used (see examples).),
#'   - or as a keyword: `"first"`, `"mid"`, `"last"` to impute to the first/mid/last
#'   day/month.
#'
#'   The argument is ignored if `highest_imputation` is less then `"D"`.
#'
#'   *Default*: `"first"`
#'
#' @param min_dates Minimum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not
#' before any of the specified dates, e.g., that the imputed adverse event start
#' date is not before the first treatment date. Only dates which are in the
#' range of possible dates of the `dtc` value are considered. The possible dates
#' are defined by the missing parts of the `dtc` date (see example below). This
#' ensures that the non-missing parts of the `dtc` date are not changed.
#' A date or date-time object is expected.
#' For example
#'
#' ```{r echo=TRUE, eval=FALSE}
#' impute_dtc_dtm(
#'   "2020-11",
#'   min_dates = list(
#'    ymd_hms("2020-12-06T12:12:12"),
#'    ymd_hms("2020-11-11T11:11:11")
#'   ),
#'   highest_imputation = "M"
#' )
#' ```
#'
#' returns `"2020-11-11T11:11:11"` because the possible dates for `"2020-11"`
#' range from `"2020-11-01T00:00:00"` to `"2020-11-30T23:59:59"`. Therefore
#' `"2020-12-06T12:12:12"` is ignored. Returning `"2020-12-06T12:12:12"` would
#' have changed the month although it is not missing (in the `dtc` date).
#'
#' @param max_dates Maximum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not after
#' any of the specified dates, e.g., that the imputed date is not after the data
#' cut off date. Only dates which are in the range of possible dates are
#' considered. A date or date-time object is expected.
#'
#' @param preserve Preserve day if month is missing and day is present
#'
#' For example `"2019---07"` would return `"2019-06-07` if `preserve = TRUE`
#' (and `date_imputation = "MID"`).
#'
#' Permitted Values: `TRUE`, `FALSE`
#'
#' Default: `FALSE`
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return A character vector
#'
#' @author Samia Kabi, Stefan Bundfuss
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' dates <- c(
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18T15",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' # No date imputation (highest_imputation defaulted to "n")
#' impute_dtc_dt(dtc = dates)
#'
#' # Impute to first day/month if date is partial
#' impute_dtc_dt(
#'   dtc = dates,
#'   highest_imputation = "M"
#' )
#' # Same as above
#' impute_dtc_dtm(
#'   dtc = dates,
#'   highest_imputation = "M",
#'   date_imputation = "01-01"
#' )
#'
#' # Impute to last day/month if date is partial
#' impute_dtc_dt(
#'   dtc = dates,
#'   highest_imputation = "M",
#'   date_imputation = "last",
#' )
#'
#' # Impute to mid day/month if date is partial
#' impute_dtc_dt(
#'   dtc = dates,
#'   highest_imputation = "M",
#'   date_imputation = "mid"
#' )
#'
#' # Impute a date and ensure that the imputed date is not before a list of
#' # minimum dates
#' impute_dtc_dt(
#'   "2020-12",
#'   min_dates = list(
#'     ymd("2020-12-06"),
#'     ymd("2020-11-11")
#'   ),
#'   highest_imputation = "M"
#' )
#'
#' # Impute completely missing dates (only possible if min_dates or max_dates is specified)
#' impute_dtc_dt(
#'   c("2020-12", NA_character_),
#'   min_dates = list(
#'     ymd("2020-12-06"),
#'     ymd("2020-11-11")
#'   ),
#'   highest_imputation = "Y"
#' )
impute_dtc_dt <- function(dtc,
                          highest_imputation = "n",
                          date_imputation = "first",
                          min_dates = NULL,
                          max_dates = NULL,
                          preserve = FALSE) {
  # Check arguments ----
  assert_character_vector(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)
  imputation_levels <- c(
    none = "n",
    day = "D",
    month = "M",
    year = "Y"
  )
  assert_character_scalar(highest_imputation, values = imputation_levels)
  highest_imputation <- dt_level(highest_imputation)
  date_imputation <-
    assert_character_scalar(
      date_imputation,
      case_sensitive = FALSE
    )
  assert_logical_scalar(preserve)

  # Parse character date ----
  two <- "(\\d{2}|-?)"
  partialdate <- stringr::str_match(dtc, paste0(
    "(\\d{4}|-?)-?",
    two,
    "-?",
    two
  ))
  partial <- vector("list", 3)
  components <- c("year", "month", "day")
  names(partial) <- components
  for (i in seq_along(components)) {
    partial[[i]] <- partialdate[, i + 1]
    partial[[i]] <- if_else(partial[[i]] %in% c("-", ""), NA_character_, partial[[i]])
  }

  # Handle preserve argument ----
  if (!preserve) {
    for (i in 2:3) {
      partial[[i]] <- if_else(is.na(partial[[i - 1]]), NA_character_, partial[[i]])
    }
  }
  # Determine target components ----
  target <- get_imputation_target_date(
    date_imputation = date_imputation,
    month = partial[["month"]]
  )

  for (c in components) {
    if (highest_imputation < dt_level(imputation_levels[[c]])) {
      target[[c]] <- "xx"
    }
  }

  # Impute ----
  imputed <- vector("list", 3)
  names(imputed) <- components
  for (c in components) {
    imputed[[c]] <- if_else(is.na(partial[[c]]), target[[c]], partial[[c]])
  }

  imputed_dtc <-
    paste(imputed[["year"]], imputed[["month"]], imputed[["day"]], sep = "-")

  imputed_dtc <-
    if_else(
      str_detect(imputed_dtc, "x"),
      NA_character_,
      imputed_dtc
    )

  if (date_imputation == "last") {
    imputed_dtc <-
      if_else(
        is.na(partial[["day"]]),
        strftime(
          rollback(ymd(imputed_dtc) + months(1)),
          format = "%Y-%m-%d"
        ),
        imputed_dtc
      )
  }

  # Handle min_dates and max_dates argument ----
  restrict_imputed_dtc_dt(
    dtc,
    imputed_dtc = imputed_dtc,
    min_dates = min_dates,
    max_dates = max_dates
  )
}

#' Create a `dt_level` object
#'
#' @param level Date level
#'
#'   *Permitted Values*: `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"n"` (none, lowest level)
#'
#' @returns A `dt_level` object
#'
#' @details A `dt_level` object is an ordered factor, i.e., two objects can be
#'   compared.
#'
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#' @keywords utils_impute
#'
dt_level <- function(level) {
  out <-
    factor(
      level,
      levels = c("n", "D", "M", "Y"),
      ordered = TRUE
    )
  class(out) <- c("dt_level", class(out))
  out
}

#' Restrict Imputed DTC date to Minimum/Maximum Dates
#'
#' @param imputed_dtc The imputed DTC date
#'
#' @inheritParams impute_dtc_dt
#'
#' @returns
#'   - The last of the minimum dates (`min_dates`) which are in the range of the
#'   partial DTC date (`dtc`)
#'   - The first of the maximum dates (`max_dates`) which are in the range of the
#'   partial DTC date (`dtc`)
#'   - `imputed_dtc` if the partial DTC date (`dtc`) is not in range of any of
#'   the minimum or maximum dates.
#'
#' @author Stefan Bundfuss
#'
#' @family utils_impute
#'
#' @keywords utils_impute
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
restrict_imputed_dtc_dt <- function(dtc,
                                    imputed_dtc,
                                    min_dates,
                                    max_dates) {
  if (!is.null(min_dates) | !is.null(max_dates)) {
    # determine range of possible dates
    min_dtc <-
      impute_dtc_dt(
        dtc,
        highest_imputation = "Y",
        date_imputation = "first"
      )
    max_dtc <-
      impute_dtc_dt(
        dtc,
        highest_imputation = "Y",
        date_imputation = "last"
      )
  }
  if (!is.null(min_dates)) {
    # for each minimum date within the range ensure that the imputed date is not
    # before it
    for (min_date in min_dates) {
      assert_date_vector(min_date)
      min_date_iso <- strftime(min_date, format = "%Y-%m-%d", tz = "UTC")
      imputed_dtc <- if_else(
        min_dtc <= min_date_iso & min_date_iso <= max_dtc,
        pmax(imputed_dtc, min_date_iso),
        imputed_dtc,
        missing = imputed_dtc
      )
    }
  }
  if (!is.null(max_dates)) {
    # for each maximum date within the range ensure that the imputed date is not
    # after it
    for (max_date in max_dates) {
      assert_date_vector(max_date)
      max_date_iso <- strftime(max_date, format = "%Y-%m-%d", tz = "UTC")
      imputed_dtc <- if_else(
        min_dtc <= max_date_iso & max_date_iso <= max_dtc,
        pmin(imputed_dtc, max_date_iso),
        imputed_dtc,
        missing = imputed_dtc
      )
    }
  }
  imputed_dtc
}

#' Convert a Date Character Vector into a Date Object
#'
#' Convert a date character vector (usually '--DTC') into a Date vector (usually '--DT').
#'
#' @param dtc The --DTC date to convert.
#'
#' @inheritParams impute_dtc_dt
#'
#' @author Samia Kabi
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return a date object
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' convert_dtc_to_dt("2019-07-18")
#' convert_dtc_to_dt("2019-07")
convert_dtc_to_dt <- function(dtc,
                              highest_imputation = "n",
                              date_imputation = "first",
                              min_dates = NULL,
                              max_dates = NULL,
                              preserve = FALSE) {
  assert_character_vector(dtc)
  warn_if_invalid_dtc(dtc, is_valid_dtc(dtc))

  imputed_dtc <- impute_dtc_dt(
    dtc = dtc,
    highest_imputation = highest_imputation,
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates,
    preserve = preserve
  )
  ymd(imputed_dtc)
}

#' Convert a Date Character Vector into a Datetime Object
#'
#' Convert a date character vector (usually `'--DTC'`) into a Date vector (usually `'--DTM'`).
#'
#' @param dtc The `'--DTC'` date to convert.
#'
#' @inheritParams impute_dtc_dtm
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return A datetime object
#'
#' @author Samia Kabi, Stefan Bundfuss
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' convert_dtc_to_dtm("2019-07-18T15:25:00")
#' convert_dtc_to_dtm("2019-07-18T00:00:00") # note Time = 00:00:00 is not printed
#' convert_dtc_to_dtm("2019-07-18")
convert_dtc_to_dtm <- function(dtc,
                               highest_imputation = "h",
                               date_imputation = "first",
                               time_imputation = "first",
                               min_dates = NULL,
                               max_dates = NULL,
                               preserve = FALSE) {
  assert_character_vector(dtc)
  warn_if_invalid_dtc(dtc, is_valid_dtc(dtc))

  dtc %>%
    impute_dtc_dtm(
      highest_imputation = highest_imputation,
      date_imputation = date_imputation,
      time_imputation = time_imputation,
      min_dates = min_dates,
      max_dates = max_dates,
      preserve = preserve
    ) %>%
    ymd_hms()
}

#' Convert a Date into a Datetime Object
#'
#' Convert a date (datetime, date, or date character) into a Date vector (usually `'--DTM'`).
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
#' @author Samia Kabi
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
  if (lubridate::is.POSIXct(dt)) {
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

#' Derive the Date Imputation Flag
#'
#' Derive the date imputation flag (`'--DTF'`) comparing a date character vector
#' (`'--DTC'`) with a Date vector (`'--DT'`).
#'
#' @param dtc The date character vector (`'--DTC'`).
#'
#'   A character date is expected in a format like `yyyy-mm-ddThh:mm:ss` (partial or complete).
#'
#' @param dt The  Date vector to compare.
#'
#'   A date object is expected.
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The date imputation flag (`'--DTF'`) (character value of `'D'`, `'M'` , `'Y'` or `NA`)
#'
#' @author Samia Kabi
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' compute_dtf(dtc = "2019-07", dt = as.Date("2019-07-18"))
#' compute_dtf(dtc = "2019", dt = as.Date("2019-07-18"))
compute_dtf <- function(dtc, dt) {
  assert_character_vector(dtc)
  assert_date_vector(dt)

  is_na <- is.na(dt)
  n_chr <- nchar(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  case_when(
    (!is_na & n_chr >= 10 & valid_dtc) | is_na | !valid_dtc ~ NA_character_,
    n_chr < 4 ~ "Y",
    n_chr == 4 ~ "M",
    n_chr == 7 ~ "D",
    n_chr == 9 ~ "M" # dates like "2019---07"
  )
}

#' Derive the Time Imputation Flag
#'
#' Derive the time imputation flag (`'--TMF'`) comparing a date character vector
#' (`'--DTC'`) with a Datetime vector (`'--DTM'`).
#'
#' @param dtc The date character vector (`'--DTC'`).
#'
#'   A character date is expected in a format like `yyyy-mm-ddThh:mm:ss` (partial or complete).
#'
#' @param dtm The Date vector to compare (`'--DTM'`).
#'
#'   A datetime object is expected.
#'
#' @param ignore_seconds_flag  ADaM IG states that given SDTM (`'--DTC'`) variable,
#' if only hours and minutes are ever collected, and seconds are imputed in
#' (`'--DTM'`) as 00, then it is not necessary to set (`'--TMF'`) to `'S'`. A user can set this
#' to `TRUE` so the `'S'` Flag is dropped from (`'--TMF'`).
#'
#'  A logical value
#'
#'   Default: `FALSE`
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The time imputation flag (`'--TMF'`) (character value of `'H'`, `'M'` , `'S'` or `NA`)
#'
#' @author Samia Kabi, Stefan Bundfuss
#'
#' @family com_date_time
#'
#' @keywords com_date_time
#'
#' @export
#'
#' @examples
#' compute_tmf(dtc = "2019-07-18T15:25", dtm = as.POSIXct("2019-07-18T15:25:00"))
#' compute_tmf(dtc = "2019-07-18T15", dtm = as.POSIXct("2019-07-18T15:25:00"))
#' compute_tmf(dtc = "2019-07-18", dtm = as.POSIXct("2019-07-18"))
compute_tmf <- function(dtc,
                        dtm,
                        ignore_seconds_flag = FALSE) {
  assert_date_vector(dtm)
  assert_character_vector(dtc)
  assert_logical_scalar(ignore_seconds_flag)

  partial <- get_partialdatetime(dtc)
  highest_miss <- convert_blanks_to_na(vector("character", length(dtc)))
  for (c in c("hour", "minute", "second")) {
    highest_miss <-
      if_else(is.na(partial[[c]]) & is.na(highest_miss), c, highest_miss)
  }
  is_na <- is.na(dtm)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  map <- c(
    hour = "H",
    minute = "M",
    second = "S"
  )
  flag <- if_else(is.na(dtm) | is.na(highest_miss), NA_character_, map[highest_miss])

  if (ignore_seconds_flag) {
    if (any(!is.na(partial[["second"]]))) {
      abort("Seconds detected in data while ignore_seconds_flag is invoked")
    } else {
      flag <- if_else(flag == "S", NA_character_, flag)
    }
  }
  flag
}

#' Derive/Impute a Date from a Date Character Vector
#'
#' Derive a date (`'--DT'`) from a date character vector (`'--DTC`').
#' The date can be imputed (see `date_imputation` argument)
#' and the date imputation flag ('`--DTF'`) can be added.
#'
#' In {admiral} we don't allow users to pick any single part of the date/time to
#' impute, we only enable to impute up to a highest level, i.e. you couldn't
#' choose to say impute months, but not days.
#'
#' @param dataset Input dataset.
#'
#'   The date character vector (`dtc`) must be present.
#'
#' @param new_vars_prefix Prefix used for the output variable(s).
#'
#'   A character scalar is expected. For the date variable "DT" is appended to
#'   the specified prefix and for the date imputation flag "DTF". I.e., for
#'   `new_vars_prefix = "AST"` the variables `ASTDT` and `ASTDTF` are created.
#'
#' @param flag_imputation Whether the date imputation flag must also be derived.
#'
#'   If `"auto"` is specified, the date imputation flag is derived if the
#'   `date_imputation` argument is not null.
#'
#'   *Default*: `"auto"`
#'
#'   *Permitted Values*: `"auto"`, `"date"` or `"none"`
#'
#'
#' @inheritParams impute_dtc_dt
#'
#' @return
#' The input dataset with the date `'--DT'` (and the date imputation flag `'--DTF'`
#' if requested) added.
#'
#' @details
#' The presence of a `'--DTF'` variable is checked and if it already exists in the input dataset,
#' a warning is issued and `'--DTF'` will be overwritten.
#'
#' @author Samia Kabi
#'
#' @family der_date_time
#'
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#'
#' mhdt <- tribble(
#'   ~MHSTDTC,
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' # Create ASTDT and ASTDTF
#' # No imputation for partial date
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC
#' )
#'
#' # Create ASTDT and ASTDTF
#' # Impute partial dates to first day/month
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M"
#' )
#'
#' # Impute partial dates to 6th of April
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "04-06"
#' )
#'
#' # Create AENDT and AENDTF
#' # Impute partial dates to last day/month
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AEN",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "last"
#' )
#'
#' # Create BIRTHDT
#' # Impute partial dates to 15th of June. No DTF
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "BIRTH",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   flag_imputation = "none"
#' )
#'
#' # Impute AE start date to the first date and ensure that the imputed date
#' # is not before the treatment start date
#' adae <- tribble(
#'   ~AESTDTC, ~TRTSDTM,
#'   "2020-12", ymd_hms("2020-12-06T12:12:12"),
#'   "2020-11", ymd_hms("2020-12-06T12:12:12")
#' )
#'
#' derive_vars_dt(
#'   adae,
#'   dtc = AESTDTC,
#'   new_vars_prefix = "AST",
#'   highest_imputation = "M",
#'   min_dates = vars(TRTSDTM)
#' )
#'
#' # A user imputing dates as middle month/day, i.e. date_imputation = "mid" can
#' # use preserve argument to "preserve" partial dates.  For example, "2019---07",
#' # will be displayed as "2019-06-07" rather than 2019-06-15 with preserve = TRUE
#'
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   preserve = TRUE
#' )
derive_vars_dt <- function(dataset,
                           new_vars_prefix,
                           dtc,
                           highest_imputation = "n",
                           date_imputation = "first",
                           flag_imputation = "auto",
                           min_dates = NULL,
                           max_dates = NULL,
                           preserve = FALSE) {
  # check and quote arguments
  assert_character_scalar(new_vars_prefix)
  assert_vars(max_dates, optional = TRUE)
  assert_vars(min_dates, optional = TRUE)
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, required_vars = vars(!!dtc))
  assert_character_scalar(
    flag_imputation,
    values = c("auto", "date", "none"),
    case_sensitive = FALSE
  )

  # output varname
  dt <- paste0(new_vars_prefix, "DT")
  warn_if_vars_exist(dataset, dt)

  # derive --DT var
  dataset <- dataset %>%
    mutate(
      !!sym(dt) := convert_dtc_to_dt(
        dtc = !!dtc,
        highest_imputation = highest_imputation,
        date_imputation = date_imputation,
        min_dates = lapply(min_dates, eval_tidy, data = rlang::as_data_mask(.)),
        max_dates = lapply(max_dates, eval_tidy, data = rlang::as_data_mask(.)),
        preserve = preserve
      )
    )

  # derive DTF
  if (flag_imputation == "date" ||
    flag_imputation == "auto" && highest_imputation != "n") {
    # add --DTF if not there already
    dtf <- paste0(new_vars_prefix, "DTF")
    dtf_exist <- dtf %in% colnames(dataset)
    if (!dtf_exist) {
      dataset <- dataset %>%
        mutate(!!sym(dtf) := compute_dtf(dtc = !!dtc, dt = !!sym(dt)))
    } else {
      msg <- sprintf(
        "The %s variable is already present in the input dataset and will not be re-derived.",
        dtf
      )
      inform(msg)
    }
  }

  dataset
}

#' Derive/Impute a Datetime from a Date Character Vector
#'
#' Derive a datetime object (`'--DTM'`) from a date character vector (`'--DTC'`).
#' The date and time can be imputed (see `date_imputation`/`time_imputation` arguments)
#' and the date/time imputation flag (`'--DTF'`, `'--TMF'`) can be added.
#'
#' In {admiral} we don't allow users to pick any single part of the date/time to
#' impute, we only enable to impute up to a highest level, i.e. you couldn't
#' choose to say impute months, but not days.
#'
#' @param dataset Input dataset
#'
#'   The date character vector (`dtc`) must be present.
#'
#' @param new_vars_prefix Prefix used for the output variable(s).
#'
#'   A character scalar is expected. For the date variable "DT" is appended to
#'   the specified prefix, for the date imputation flag "DTF", and for the time
#'   imputation flag "TMF". I.e., for `new_vars_prefix = "AST"` the variables
#'   `ASTDT`, `ASTDTF`, and `ASTTMF` are created.
#'
#'
#' @param flag_imputation Whether the date/time imputation flag(s) must also be derived.
#'
#'   If `"auto"` is specified, the date imputation flag is derived if the
#'   `date_imputation` argument is not null and the time imputation flag is
#'   derived if the `time_imputation` argument is not null
#'
#'   *Default*: `"auto"`
#'
#'   *Permitted Values*: `"auto"`, `"date"`, `"time"`, `"both"`, or `"none"`
#'
#'
#' @inheritParams impute_dtc_dtm
#' @inheritParams compute_tmf
#'
#' @details
#' The presence of a `'--DTF'` variable is checked and the variable is not derived
#' if it already exists in the input dataset. However, if `'--TMF'` already exists
#' in the input dataset, a warning is issued and `'--TMF'` will be overwritten.
#'
#' @return  The input dataset with the datetime `'--DTM'` (and the date/time imputation
#' flag `'--DTF'`, `'--TMF'`) added.
#'
#' @author Samia Kabi
#'
#' @family der_date_time
#'
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#'
#' mhdt <- tribble(
#'   ~MHSTDTC,
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M"
#' )
#'
#' # Impute AE end date to the last date and ensure that the imputed date is not
#' # after the death or data cut off date
#' adae <- tribble(
#'   ~AEENDTC, ~DTHDT, ~DCUTDT,
#'   "2020-12", ymd("2020-12-06"), ymd("2020-12-24"),
#'   "2020-11", ymd("2020-12-06"), ymd("2020-12-24")
#' )
#'
#' derive_vars_dtm(
#'   adae,
#'   dtc = AEENDTC,
#'   new_vars_prefix = "AEN",
#'   highest_imputation = "M",
#'   date_imputation = "last",
#'   time_imputation = "last",
#'   max_dates = vars(DTHDT, DCUTDT)
#' )
#'
#' # Seconds has been removed from the input dataset.  Function now uses
#' # ignore_seconds_flag to remove the 'S' from the --TMF variable.
#' mhdt <- tribble(
#'   ~MHSTDTC,
#'   "2019-07-18T15:25",
#'   "2019-07-18T15:25",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   ignore_seconds_flag = TRUE
#' )
#'
#' # A user imputing dates as middle month/day, i.e. date_imputation = "MID" can
#' # use preserve argument to "preserve" partial dates.  For example, "2019---07",
#' # will be displayed as "2019-06-07" rather than 2019-06-15 with preserve = TRUE
#'
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   preserve = TRUE
#' )
derive_vars_dtm <- function(dataset,
                            new_vars_prefix,
                            dtc,
                            highest_imputation = "h",
                            date_imputation = "first",
                            time_imputation = "first",
                            flag_imputation = "auto",
                            min_dates = NULL,
                            max_dates = NULL,
                            preserve = FALSE,
                            ignore_seconds_flag = FALSE) {
  # check and quote arguments
  assert_character_scalar(new_vars_prefix)
  assert_vars(max_dates, optional = TRUE)
  assert_vars(min_dates, optional = TRUE)
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, required_vars = vars(!!dtc))
  assert_character_scalar(
    flag_imputation,
    values = c("auto", "both", "date", "time", "none"),
    case_sensitive = FALSE
  )

  dtm <- paste0(new_vars_prefix, "DTM")

  # Issue a warning if --DTM already exists
  warn_if_vars_exist(dataset, dtm)
  mask <- rlang::as_data_mask(dataset)
  dataset[[dtm]] <- convert_dtc_to_dtm(
    dtc = eval_tidy(dtc, mask),
    highest_imputation = highest_imputation,
    date_imputation = date_imputation,
    time_imputation = time_imputation,
    min_dates = lapply(min_dates, eval_tidy, data = mask),
    max_dates = lapply(max_dates, eval_tidy, data = mask),
    preserve = preserve
  )

  if (flag_imputation %in% c("both", "date") ||
    flag_imputation == "auto" && dtm_level(highest_imputation) > dtm_level("h")) {
    # add --DTF if not there already
    dtf <- paste0(new_vars_prefix, "DTF")
    dtf_exist <- dtf %in% colnames(dataset)
    if (!dtf_exist) {
      dataset <- dataset %>%
        mutate(!!sym(dtf) := compute_dtf(dtc = !!dtc, dt = !!sym(dtm)))
    } else {
      msg <- sprintf(
        "The %s variable is already present in the input dataset and will not be re-derived.",
        dtf
      )
      inform(msg)
    }
  }

  if (flag_imputation %in% c("both", "time") ||
    flag_imputation == "auto" && highest_imputation != "n") {
    # add --TMF variable
    tmf <- paste0(new_vars_prefix, "TMF")
    warn_if_vars_exist(dataset, tmf)

    dataset <- dataset %>%
      mutate(!!sym(tmf) := compute_tmf(
        dtc = !!dtc,
        dtm = !!sym(dtm),
        ignore_seconds_flag = ignore_seconds_flag
      ))
  }

  dataset
}
