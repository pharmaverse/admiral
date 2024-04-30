#' Derive/Impute a Datetime from a Date Character Vector
#'
#' Derive a datetime object (`'--DTM'`) from a date character vector (`'--DTC'`).
#' The date and time can be imputed (see `date_imputation`/`time_imputation` arguments)
#' and the date/time imputation flag (`'--DTF'`, `'--TMF'`) can be added.
#'
#' In `{admiral}` we don't allow users to pick any single part of the date/time to
#' impute, we only enable to impute up to a highest level, i.e. you couldn't
#' choose to say impute months, but not days.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("dtc"))`
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
#'   If `"both"` or `"date"` is specified, then date imputation flag is derived.
#'   If `"auto"` is specified and `highest_imputation` argument is greater than
#'   `"h"`, then date imputation flag is derived.
#'
#'   If `"both"` or `"time"` is specified, then time imputation flag is derived.
#'   If `"auto"` is specified and `highest_imputation` argument is not `"n"`,
#'   then time imputation flag is derived.
#'
#'   If `"none"` is specified, then no date or time imputation flag is derived.
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
#'   max_dates = exprs(DTHDT, DCUTDT)
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
derive_vars_dtm <- function(dataset, # nolint: cyclocomp_linter
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
  dtc <- assert_symbol(enexpr(dtc))
  assert_data_frame(dataset, required_vars = exprs(!!dtc))
  assert_character_scalar(
    flag_imputation,
    values = c("auto", "both", "date", "time", "none"),
    case_sensitive = FALSE
  )
  if ((highest_imputation == "Y" && is.null(min_dates) && is.null(max_dates)) ||
    (highest_imputation == "Y" && length(min_dates) == 0 && length(max_dates) == 0)) {
    cli_abort(paste(
      "If {.code highest_impuation = \"Y\"} is specified, {.arg min_dates} or",
      "{.arg max_dates} must be specified respectively."
    ))
  }
  if (highest_imputation == "Y") {
    assert_character_scalar(date_imputation, values = c("first", "last"))
  }
  if (highest_imputation == "Y" && is.null(min_dates) && date_imputation == "first") {
    cli_warn(paste(
      "If {.code highest_impuation = \"Y\"} and {.code date_imputation = \"first\"}",
      "is specified, {.arg min_dates} should be specified."
    ))
  }
  if (highest_imputation == "Y" && is.null(max_dates) && date_imputation == "last") {
    cli_warn(paste(
      "If {.code highest_impuation = \"Y\"} and {.code date_imputation = \"last\"}",
      "is specified, {.arg max_dates} should be specified."
    ))
  }

  dtm <- paste0(new_vars_prefix, "DTM")

  # Issue a warning if --DTM already exists
  warn_if_vars_exist(dataset, dtm)
  mask <- as_data_mask(dataset)

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
      cli_inform(paste(
        "The {.var {dtf}} variable is already present in the input dataset and",
        "will not be re-derived."
      ))
    }
  }

  if (flag_imputation %in% c("both", "time") ||
    flag_imputation == "auto" && highest_imputation != "n") {
    # add --TMF variable
    tmf <- paste0(new_vars_prefix, "TMF")
    warn_if_vars_exist(dataset, tmf)

    tryCatch(
      dataset <- dataset %>%
        mutate(!!sym(tmf) := compute_tmf(
          dtc = !!dtc,
          dtm = !!sym(dtm),
          ignore_seconds_flag = ignore_seconds_flag
        )),
      "dplyr:::mutate_error" = function(cnd) {
        cli_abort(
          message = cnd$parent$message,
          call = parent.frame(n = 4)
        )
      }
    )
  }

  dataset
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

  imputed_dtc <- impute_dtc_dtm(
    dtc = dtc,
    highest_imputation = highest_imputation,
    date_imputation = date_imputation,
    time_imputation = time_imputation,
    min_dates = min_dates,
    max_dates = max_dates,
    preserve = preserve
  )

  imputed_dtc <- if_else(
    str_starts(imputed_dtc, "(0000|9999)") | imputed_dtc %in% c("0000-01-01", "9999-12-31"), # nolint
    NA_character_,
    imputed_dtc
  )

  ymd_hms(imputed_dtc)
}

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
#'   *Permitted Values*: `"Y"` (year, highest level), `"M"` (month), `"D"`
#'   (day), `"h"` (hour), `"m"` (minute), `"s"` (second), `"n"` (none, lowest
#'   level)
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

#' @param preserve Preserve lower level date/time part when higher order part
#' is missing, e.g. preserve day if month is missing or
#' preserve minute when hour is missing.
#'
#' For example `"2019---07"` would return `"2019-06-07` if `preserve = TRUE`
#' (and `date_imputation = "mid"`).
#'
#' Permitted Values: `TRUE`, `FALSE`
#'
#' @inheritParams impute_dtc_dt
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return A character vector
#'
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
#'     ymd_hms("2020-12-06T12:12:12", "2020-01-01T01:01:01"),
#'     ymd_hms("2020-11-11T11:11:11", NA)
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
  restricted <- restrict_imputed_dtc_dtm(
    dtc,
    imputed_dtc = imputed_dtc,
    min_dates = min_dates,
    max_dates = max_dates
  )

  if (highest_imputation == "Y" && is.null(min_dates) && is.null(max_dates)) {
    warning("If `highest_impuation` = \"Y\" is specified, `min_dates` or `max_dates` should be specified respectively.") # nolint
  }

  return(restricted)
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
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
restrict_imputed_dtc_dtm <- function(dtc,
                                     imputed_dtc,
                                     min_dates,
                                     max_dates) {
  if (!(is.null(min_dates) || length(min_dates) == 0) ||
    !(is.null(max_dates) || length(max_dates) == 0)) {
    suppress_warning(
      { # nolint
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
      },
      # Suppress warning because we need to run without min/max dates but users should not
      regexpr = "If `highest_impuation` = \"Y\" is specified, `min_dates` or `max_dates` should be specified respectively." # nolint
    )
  }
  if (!(is.null(min_dates) || length(min_dates) == 0)) {
    if (length(unique(c(length(imputed_dtc), unlist(lapply(min_dates, length))))) != 1) {
      cli_abort("Length of {.arg min_dates} do not match length of dates to be imputed.")
    }
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
  if (!(is.null(max_dates) || length(max_dates) == 0)) {
    if (length(unique(c(length(imputed_dtc), unlist(lapply(max_dates, length))))) != 1) {
      cli_abort("Length of {.arg max_dates} do not match length of dates to be imputed.")
    }
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
#' *Permitted Values*: A logical value
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The time imputation flag (`'--TMF'`) (character value of `'H'`, `'M'` , `'S'` or `NA`)
#'
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
#' compute_tmf(dtc = "2019-07-18T15:25", dtm = ymd_hms("2019-07-18T15:25:00"))
#' compute_tmf(dtc = "2019-07-18T15", dtm = ymd_hms("2019-07-18T15:25:00"))
#' compute_tmf(dtc = "2019-07-18", dtm = ymd("2019-07-18"))
#' compute_tmf(dtc = "2022-05--T00:00", dtm = ymd_hms("2022-05-15T23:59:59"))
#' compute_tmf(dtc = "2022-05--T23:00", dtm = ymd_hms("2022-05-15T23:59:59"))
#' compute_tmf(dtc = "2022-05--T23:59:00", dtm = ymd_hms("2022-05-15T23:59:59"))
#'
compute_tmf <- function(dtc,
                        dtm,
                        ignore_seconds_flag = FALSE) {
  assert_date_vector(dtm)
  assert_character_vector(dtc)
  assert_logical_scalar(ignore_seconds_flag)

  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  partial <- get_partialdatetime(dtc)
  highest_miss <- convert_blanks_to_na(vector("character", length(dtc)))

  # concatenate lubridate functions: `hour()`, `minute()`, `second()` to map over dtm input
  hms <- c("hour", "minute", "second")

  # extract hour, minute, second over each value of dtm and put into a list time_part
  time_part <-
    map(set_names(hms), function(y) map_dbl(dtm, function(x) exec(y, x)))

  for (c in hms) {
    highest_miss <-
      if_else((is.na(partial[[c]]) & is.na(highest_miss)) |
        (
          !is.na(partial[[c]]) &
            is.na(highest_miss) & as.numeric(partial[[c]]) != time_part[[c]]
        ),
      c,
      highest_miss
      )
  }

  map <- c(
    hour = "H",
    minute = "M",
    second = "S"
  )
  flag <- if_else(is.na(dtm) | is.na(highest_miss), NA_character_, unname(map[highest_miss]))

  if (ignore_seconds_flag) {
    if (any(!is.na(partial[["second"]]))) {
      cli_abort(
        "Seconds detected in data while {.arg ignore_seconds_flag} is invoked"
      )
    } else {
      flag <- if_else(flag == "S", NA_character_, flag)
    }
  }
  flag
}
