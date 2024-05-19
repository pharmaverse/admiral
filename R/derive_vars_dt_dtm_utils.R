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
#'
#' @family utils_impute
#' @keywords internal
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

  return(target)
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

  return(target)
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
