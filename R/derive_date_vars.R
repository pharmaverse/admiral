#' Impute Partial Date(-time) Portion of a `'--DTC'` Variable
#'
#' Imputation partial date/time portion of a `'--DTC'` variable. based on user
#' input.
#'
#' @param dtc The `'--DTC'` date to impute
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. If the year part is not recorded (missing date), no
#'   imputation is performed.
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   If `NULL`: no date imputation is performed and partial dates are returned as
#'   missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for the 15th
#'   of June,
#'   - or as a keyword: `"FIRST"`, `"MID"`, `"LAST"` to impute to the first/mid/last
#'   day/month.
#'
#'   Default is `NULL`.
#'
#' @param time_imputation The value to impute the time when a timepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as `"hh:mm:ss"`: e.g. `"00:00:00"`
#'   for the start of the day,
#'   - or as a keyword: `"FIRST"`,`"LAST"` to impute to the start/end of a day.
#'
#'   Default is `"00:00:00"`.
#'
#'
#'
#'
#' @param min_dates Minimum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not
#' before any of the specified dates, e.g., that the imputed adverse event start
#' date is not before the first treatment date. Only dates which are in the
#' range of possible dates of the `dtc` value are considered. The possible dates
#' are defined by the missing parts of the `dtc` date (see example below). This
#' ensures that the non-missing parts of the `dtc` date are not changed. For
#' example
#'
#' ```
#' impute_dtc(
#'   "2020-11",
#'   min_dates = list(
#'     ymd_hms("2020-12-06T12:12:12"),
#'     ymd_hms("2020-11-11T11:11:11")
#'    ),

#'   date_imputation = "first"
#' )
#' ```
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
#' considered.
#'
#' @param preserve Preserve partial dates when doing date imputation for middle
#' day and month
#'
#' A user wishing to preserve partial dates when doing middle day and month date
#' imputation can invoke this argument.  For example `"2019---07"` would return
#' `"2019-06-07` if date_imputation = "MID" and preserve = TRUE.
#'
#'  A logical value
#'
#'  Default: `FALSE`
#'
#'
#' @author Samia Kabi
#'
#' @return A character vector
#'
#' @keywords computation timing
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
#' # No date imputation (date_imputation defaulted to NULL)
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc(dtc = dates)
#'
#' # No date imputation (date_imputation defaulted to NULL)
#' # Missing time part imputed with 23:59:59 portion
#' impute_dtc(
#'   dtc = dates,
#'   time_imputation = "23:59:59"
#' )
#'
#' # Same as above
#' impute_dtc(
#'   dtc = dates,
#'   time_imputation = "LAST"
#' )
#'
#' # Impute to first day/month if date is partial
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc(
#'   dtc = dates,
#'   date_imputation = "01-01"
#' )
#' # same as above
#' impute_dtc(
#'   dtc = dates,
#'   date_imputation = "FIRST"
#' )
#'
#' # Impute to last day/month if date is partial
#' # Missing time part imputed with 23:59:59 portion
#' impute_dtc(
#'   dtc = dates,
#'   date_imputation = "LAST",
#'   time_imputation = "LAST"
#' )
#'
#' # Impute to mid day/month if date is partial
#' # Missing time part imputed with 00:00:00 portion by default
#' impute_dtc(
#'   dtc = dates,
#'   date_imputation = "MID"
#' )
#'
#' # Impute a date and ensure that the imputed date is not before a list of
#' # minimum dates
#' impute_dtc(
#'   "2020-12",
#'   min_dates = list(
#'     ymd_hms("2020-12-06T12:12:12"),
#'     ymd_hms("2020-11-11T11:11:11")
#'   ),
#'   date_imputation = "first"
#' )
impute_dtc <- function(dtc,
                       date_imputation = NULL,
                       time_imputation = "00:00:00",
                       min_dates = NULL,
                       max_dates = NULL,
                       preserve = FALSE) {
  # Issue a warning if incorrect  DTC is present
  n_chr <- nchar(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)
  assert_logical_scalar(preserve)

  # date imputation
  if (!is.null(date_imputation)) {
    # check input for date_imputation
    assert_that(is_valid_date_entry(date_imputation))

    # Specific setup for FIRST/MID/LAST
    # make keywords case-insensitive
    date_imputation <- str_to_upper(date_imputation)
    if (date_imputation %in% c("FIRST", "MID", "LAST")) {
      d <- list(FIRST = "01", MID = "15", LAST = "01")[[date_imputation]]
      mo <- list(FIRST = "01", MID = "06", LAST = "12")[[date_imputation]]
    } else {
      # otherwise, use time_imputation input
      day__ <- as.integer(sub(".*-", "", date_imputation))
      mo__ <- as.integer(sub("-.*", "", date_imputation))
      # check input for day and moth are valid
      assert_that(is_valid_day(day__))
      assert_that(is_valid_month(mo__))

      d <- sprintf("%02d", day__)
      mo <- sprintf("%02d", mo__)
    }

    imputed_date <- case_when(
      !valid_dtc ~ NA_character_,
      n_chr >= 10 ~ substr(dtc, 1, 10),
      n_chr == 9 ~ paste0(substr(dtc, 1, 4), "-", mo, "-", d),
      n_chr == 7 ~ paste0(dtc, "-", d),
      n_chr == 4 ~ paste0(dtc, "-", mo, "-", d)
    )

    # 3 blocks of  if/else statements that deal with date imputation and
    # preserving partial dates.
    # Ex: 2019---07 with MID and preserve = TRUE gives 2019-06-07
    if (date_imputation == "MID" & preserve) {

      imputed_date <- case_when(
        n_chr == 9 ~ paste0(substr(dtc, 1, 4), "-", "06", "-", substr(dtc, 8, 9)),
        n_chr == 4  ~ paste0(dtc, "-", "06", "-", "30"),
        TRUE ~ imputed_date)

    } else if (date_imputation == "MID" & !preserve) {

      imputed_date <- case_when(
        n_chr == 9  ~ paste0(substr(dtc, 1, 4), "-", mo, "-", d),
        n_chr == 4  ~ paste0(dtc, "-", "06", "-", "30"),
        TRUE ~ imputed_date)

    } else if (date_imputation != "MID" & preserve) {

      imputed_date <- case_when(
        n_chr == 9  ~ paste0(substr(dtc, 1, 4), "-", mo, "-", substr(dtc, 8, 9)),
        TRUE ~ imputed_date)
    }

    # Ex: 2019---07 with LAST and preserve = TRUE gives 2019-12-07
    if (date_imputation == "LAST" & !preserve) {

      imputed_date <- case_when(
        n_chr < 10 & date_imputation == "LAST" & !preserve ~
          as.character(
            ceiling_date(
              as.Date(imputed_date, format = "%Y-%m-%d"), "month") - days(1)),
        TRUE ~ imputed_date)


    } else if (date_imputation == "LAST" & preserve) {

      imputed_date <- case_when(
        n_chr == 9 ~  paste0(substr(dtc, 1, 4), "-", "12", "-", substr(dtc, 8, 9)),
        n_chr %in% c(4, 7) ~
          as.character(
            ceiling_date(
              as.Date(imputed_date, format = "%Y-%m-%d"), "month") - days(1)),
        TRUE ~ imputed_date)
    }

    # Ex: 2019---07 with FIRST and preserve = TRUE gives 2019-01-07
    if (date_imputation == "FIRST" & preserve) {

       imputed_date <- case_when(
        n_chr == 9 & date_imputation == "FIRST" & preserve ~
          paste0(substr(dtc, 1, 4), "-", "01", "-", substr(dtc, 8, 9)),
        TRUE ~ imputed_date)
    }

  } else  {
    # no imputation
    imputed_date <- if_else(n_chr >= 10 & valid_dtc, substr(dtc, 1, 10), NA_character_)
  }

  if (!is.null(time_imputation)) {
    # impute time
    assert_that(is_valid_time_entry(time_imputation))
    # make keywords case-insensitive
    time_imputation <- str_to_upper(time_imputation)
    if (time_imputation == "FIRST") {
      imputed_time <- "00:00:00"
      sec <- ":00"
      min <- ":00"
      h <- "00"
    } else if (time_imputation == "LAST") {
      imputed_time <- "23:59:59"
      sec <- ":59"
      min <- ":59"
      h <- "23"
    } else {
      imputed_time <- time_imputation
      sec <- paste0(":", paste0(substr(dtc, 18, 19), substr(time_imputation, 7, 8)))
      min <- paste0(":", paste0(substr(dtc, 15, 16), substr(time_imputation, 4, 5)))
      h <- paste0(substr(dtc, 12, 13), substr(time_imputation, 1, 2))
    }

    imputed_time <- case_when(
      n_chr >= 19 ~ substr(dtc, 12, 19),
      n_chr == 16 ~ paste0(substr(dtc, 12, 16), sec),
      n_chr == 13 ~ paste0(substr(dtc, 12, 13), min, sec),
      n_chr == 10 ~ paste0(h, min, sec),
      TRUE ~ imputed_time
    )
  } else {
    # no imputation
    imputed_time <- if_else(n_chr >= 19 & valid_dtc, substr(dtc, 12, 19), NA_character_)
  }

  imputed_dtc <- if_else(
    !is.na(imputed_date) & !is.na(imputed_time),
    paste0(imputed_date, "T", imputed_time),
    NA_character_
  )

  # adjust imputed date to minimum and maximum dates
  if (!is.null(min_dates) | !is.null(max_dates)) {
    # determine range of possible dates
    min_dtc <- impute_dtc(dtc, date_imputation = "first", time_imputation = "first")
    max_dtc <- impute_dtc(dtc, date_imputation = "last", time_imputation = "last")
  }
  if (!is.null(min_dates)) {
    # for each minimum date within the range ensure that the imputed date is not
    # before it
    for (min_date in min_dates) {
      assert_that(is_date(min_date))
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
      assert_that(is_date(max_date))
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

#' Convert a Date Character Vector into a Date Object
#'
#' Convert a date character vector (usually '--DTC') into a Date vector (usually '--DT').
#'
#' @param dtc The --DTC date to convert.
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss.
#'   A partial date will return a NA date and a warning will be issued:
#'   'All formats failed to parse. No formats found.'.
#'   Note: you can use impute_dtc function to build a complete date.
#'
#' @inheritParams impute_dtc
#'
#' @author Samia Kabi
#'
#' @return a date object
#'
#' @keywords computation timing
#'
#' @export
#'
#' @examples
#' convert_dtc_to_dt("2019-07-18")
#' convert_dtc_to_dt("2019-07")
convert_dtc_to_dt <- function(dtc,
                              date_imputation = NULL,
                              min_dates = NULL,
                              max_dates = NULL,
                              preserve = FALSE) {
  assert_that(is.character(dtc))
  warn_if_invalid_dtc(dtc, is_valid_dtc(dtc))

  imputed_dtc <- impute_dtc(
    dtc = dtc,
    date_imputation = date_imputation,
    time_imputation = "first",
    min_dates = min_dates,
    max_dates = max_dates,
    preserve = preserve
  )

  if_else(
    is.na(imputed_dtc),
    ymd(NA),
    ymd(substr(imputed_dtc, 1, 10))
  )
}

#' Convert a Date Character Vector into a Datetime Object
#'
#' Convert a date character vector (usually `'--DTC'`) into a Date vector (usually `'--DTM'`).
#'
#' @param dtc The `'--DTC'` date to convert.
#'
#'   A character date is expected in a format like `yyyy-mm-ddThh:mm:ss`.
#'   A partial datetime will issue a warning.
#'   Note: you can use [impute_dtc()] function to build a complete datetime.
#'
#' @inheritParams impute_dtc
#'
#' @author Samia Kabi
#'
#' @return A datetime object
#'
#' @keywords computation timing
#'
#' @export
#'
#' @examples
#' convert_dtc_to_dtm("2019-07-18T15:25:00")
#' convert_dtc_to_dtm("2019-07-18T00:00:00") # note Time = 00:00:00 is not printed
#' convert_dtc_to_dtm("2019-07-18")
convert_dtc_to_dtm <- function(dtc,
                               date_imputation = NULL,
                               time_imputation = NULL,
                               min_dates = NULL,
                               max_dates = NULL,
                               preserve = FALSE) {
  assert_character_vector(dtc)
  warn_if_invalid_dtc(dtc, is_valid_dtc(dtc))

  dtc %>%
    impute_dtc(
      date_imputation = date_imputation,
      time_imputation = time_imputation,
      min_dates = min_dates,
      max_dates = max_dates,
      preserve = preserve
    ) %>%
    as_iso_dtm()
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
#' @author Samia Kabi
#'
#' @return A datetime object
#'
#' @keywords computation timing
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
                                date_imputation = NULL,
                                time_imputation = NULL,
                                min_dates = NULL,
                                max_dates = NULL,
                                preserve = FALSE) {

  if (lubridate::is.POSIXct(dt)) {
    return(dt)
  }
  else {
    if (is_date(dt)) {
      dt <- format(dt, "%Y-%m-%d")
    }

    # convert dtc to dtm
    dt %>%
      convert_dtc_to_dtm(
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
#' @author Samia Kabi
#'
#' @return The date imputation flag (`'--DTF'`) (character value of `'D'`, `'M'` , `'Y'` or `NA`)
#'
#' @keywords computation timing
#'
#' @export
#'
#' @examples
#' compute_dtf(dtc = "2019-07", dt = as.Date("2019-07-18"))
#' compute_dtf(dtc = "2019", dt = as.Date("2019-07-18"))
compute_dtf <- function(dtc, dt) {
  assert_that(is.character(dtc), is_date(dt))

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
#' @author Samia Kabi
#'
#' @return The time imputation flag (`'--TMF'`) (character value of `'H'`, `'M'` , `'S'` or `NA`)
#'
#' @keywords computation timing
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

  assert_that(is_date(dtm))
  assert_character_vector(dtc)
  assert_logical_scalar(ignore_seconds_flag)

  is_na <- is.na(dtm)
  n_chr <- nchar(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)


  if (ignore_seconds_flag)  {
    if ((any(n_chr >= 17))) {
      abort("Seconds detected in data while ignore_seconds_flag is invoked")
    } else {
      case_when(
        (!is_na & n_chr >= 19 & valid_dtc) | is_na | !valid_dtc ~ NA_character_,
        n_chr == 13 ~ "M",
        n_chr == 10 | (n_chr > 0 & n_chr < 10) ~ "H")
      }
  } else {
    case_when(
      (!is_na & n_chr >= 19 & valid_dtc) | is_na | !valid_dtc ~ NA_character_,
      n_chr == 16 ~ "S",
      n_chr == 13 ~ "M",
      n_chr == 10 | (n_chr > 0 & n_chr < 10) ~ "H")
     }
}

#' Derive/Impute a Date from a Date Character Vector
#'
#' Derive a date (`'--DT'`) from a date character vector (`'--DTC`').
#' The date can be imputed (see `date_imputation` parameter)
#' and the date imputation flag ('`--DTF'`) can be added.
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
#' @param flag_imputation Whether the date imputation flag should also be derived.
#'
#'   A logical value
#'
#'   Default: `TRUE`
#'
#'
#' @inheritParams impute_dtc
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
#' @keywords adam derivation timing
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' mhdt <- tibble::tribble(
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
#' # no imputation for partial date
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
#'   date_imputation = "FIRST"
#' )
#'
#' # Impute partial dates to 6th of April
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   date_imputation = "04-06"
#' )
#'
#' # Create AENDT and AENDTF
#' # Impute partial dates to last day/month
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AEN",
#'   dtc = MHSTDTC,
#'   date_imputation = "LAST"
#' )
#'
#' # Create BIRTHDT
#' # Impute partial dates to 15th of June. No DTF
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "BIRTH",
#'   dtc = MHSTDTC,
#'   date_imputation = "MID",
#'   flag_imputation = FALSE
#' )
#'
#' # Impute AE start date to the first date and ensure that the imputed date
#' # is not before the treatment start date
#' adae <- tibble::tribble(
#'   ~AESTDTC, ~TRTSDTM,
#'   "2020-12", ymd_hms("2020-12-06T12:12:12"),
#'   "2020-11", ymd_hms("2020-12-06T12:12:12")
#' )
#'
#' derive_vars_dt(
#'   adae,
#'   dtc = AESTDTC,
#'   new_vars_prefix = "AST",
#'   date_imputation = "first",
#'   min_dates = vars(TRTSDTM)
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
#'   date_imputation = "MID",
#'   preserve = TRUE
#' )
derive_vars_dt <- function(dataset,
                           new_vars_prefix,
                           dtc,
                           date_imputation = NULL,
                           flag_imputation = TRUE,
                           min_dates = NULL,
                           max_dates = NULL,
                           preserve = FALSE) {

  # check and quote parameters
  assert_character_scalar(new_vars_prefix)
  assert_vars(max_dates, optional = TRUE)
  assert_vars(min_dates, optional = TRUE)
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, required_vars = vars(!!dtc))
  assert_logical_scalar(flag_imputation)

  # output varname
  dt <- paste0(new_vars_prefix, "DT")
  warn_if_vars_exist(dataset, dt)

  # derive --DT var
  dataset <- dataset %>%
    mutate(
      !!sym(dt) := convert_dtc_to_dt(
        dtc = !!dtc,
        date_imputation = date_imputation,
        min_dates = lapply(min_dates, eval_tidy, data = rlang::as_data_mask(.)),
        max_dates = lapply(max_dates, eval_tidy, data = rlang::as_data_mask(.)),
        preserve = preserve
      )
    )

  # derive DTF
  if (flag_imputation) {
    dtf <- paste0(new_vars_prefix, "DTF")
    warn_if_vars_exist(dataset, dtf)
    dataset <- dataset %>%
      mutate(!!sym(dtf) := compute_dtf(dtc = !!dtc, dt = !!sym(dt)))
  }

  dataset
}

#' Derive/Impute a Datetime from a Date Character Vector
#'
#' Derive a datetime object (`'--DTM'`) from a date character vector (`'--DTC'`).
#' The date and time can be imputed (see `date_imputation`/`time_imputation` parameters)
#' and the date/time imputation flag (`'--DTF'`, `'--TMF'`) can be added.
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
#'   `date_imputation` parameter is not null and the time imputation flag is
#'   derived if the `time_imputation` parameter is not null
#'
#'   *Default*: `"auto"`
#'
#'   *Permitted Values*: `"auto"`, `"date"`, `"time"`, `"both"`, or `"none"`
#'
#'
#' @inheritParams impute_dtc
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
#' @keywords derivation adam timing
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' mhdt <- tibble::tribble(
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
#'   date_imputation = "FIRST",
#'   time_imputation = "FIRST"
#' )
#'
#' # Impute AE end date to the last date and ensure that the imputed date is not
#' # after the death or data cut off date
#' adae <- tibble::tribble(
#'   ~AEENDTC, ~DTHDT, ~DCUTDT,
#'   "2020-12", ymd("2020-12-06"), ymd("2020-12-24"),
#'   "2020-11", ymd("2020-12-06"), ymd("2020-12-24")
#' )
#'
#' derive_vars_dtm(
#'   adae,
#'   dtc = AEENDTC,
#'   new_vars_prefix = "AEN",
#'   date_imputation = "last",
#'   time_imputation = "last",
#'   max_dates = vars(DTHDT, DCUTDT)
#' )
#'
#' # Seconds has been removed from the input dataset.  Function now uses
#' # ignore_seconds_flag to remove the 'S' from the --TMF variable.
#' mhdt <- tibble::tribble(
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
#'   date_imputation = "FIRST",
#'   time_imputation = "FIRST",
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
#'   date_imputation = "MID",
#'   preserve = TRUE
#' )
derive_vars_dtm <- function(dataset,
                            new_vars_prefix,
                            dtc,
                            date_imputation = NULL,
                            time_imputation = "00:00:00",
                            flag_imputation = "auto",
                            min_dates = NULL,
                            max_dates = NULL,
                            preserve = FALSE,
                            ignore_seconds_flag = FALSE) {

  # check and quote parameters
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
    date_imputation = date_imputation,
    time_imputation = time_imputation,
    min_dates = lapply(min_dates, eval_tidy, data = mask),
    max_dates = lapply(min_dates, eval_tidy, data = mask),
    preserve = preserve
  )

  if (flag_imputation %in% c("both", "date") ||
      flag_imputation == "auto" && !is.null(date_imputation)) {
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
      flag_imputation == "auto" && !is.null(time_imputation)) {
    # add --TMF variable
    tmf <- paste0(new_vars_prefix, "TMF")
    warn_if_vars_exist(dataset, tmf)

    dataset <- dataset %>%
      mutate(!!sym(tmf) := compute_tmf(
        dtc = !!dtc,
        dtm = !!sym(dtm),
        ignore_seconds_flag = ignore_seconds_flag)
        )

  }


  dataset
}
