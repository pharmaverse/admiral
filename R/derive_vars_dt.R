#' Derive/Impute a Date from a Date Character Vector
#'
#' Derive a date (`'--DT'`) from a date character vector (`'--DTC`').
#' The date can be imputed (see `date_imputation` argument)
#' and the date imputation flag ('`--DTF'`) can be added.
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
#'   the specified prefix and for the date imputation flag "DTF". I.e., for
#'   `new_vars_prefix = "AST"` the variables `ASTDT` and `ASTDTF` are created.
#'
#' @param flag_imputation Whether the date imputation flag must also be derived.
#'
#'   If `"auto"` is specified and `highest_imputation` argument is not `"n"`,
#'   then date imputation flag is derived.
#'
#'   If `"date"` is specified, then date imputation flag is derived.
#'
#'   If `"none"` is specified, then no date imputation flag is derived.
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
#' # Impute partial dates to 15th of June. No Date Imputation Flag
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
#'   min_dates = exprs(TRTSDTM)
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
derive_vars_dt <- function(dataset, # nolint: cyclocomp_linter
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
  dtc <- assert_symbol(enexpr(dtc))
  assert_data_frame(dataset, required_vars = exprs(!!dtc))
  assert_character_scalar(
    flag_imputation,
    values = c("auto", "date", "none"),
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
        min_dates = lapply(min_dates, eval_tidy, data = as_data_mask(.)),
        max_dates = lapply(max_dates, eval_tidy, data = as_data_mask(.)),
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

#' Convert a Date Character Vector into a Date Object
#'
#' Convert a date character vector (usually '--DTC') into a Date vector (usually '--DT').
#'
#' @param dtc The --DTC date to convert.
#'
#' @inheritParams impute_dtc_dt
#'
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
  imputed_dtc <- if_else(
    str_starts(imputed_dtc, "(0000|9999)") | imputed_dtc %in% c("0000-01-01", "9999-12-31"), # nolint
    NA_character_,
    imputed_dtc
  )
  ymd(imputed_dtc)
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
#'   day/month. If `"mid"` is specified, missing components are imputed as the
#'   middle of the possible range:
#'       - If both month and day are missing, they are imputed as `"06-30"`
#'        (middle of the year).
#'       - If only day is missing, it is imputed as `"15"` (middle of the month).
#'
#'   The argument is ignored if `highest_imputation` is less then `"D"`.
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
#' # No date imputation (highest_imputation defaulted to "n")
#' impute_dtc_dt(dtc = dates)
#'
#' # Impute to first day/month if date is partial
#' impute_dtc_dt(
#'   dtc = dates,
#'   highest_imputation = "M"
#' )
#' # Same as above
#' impute_dtc_dt(
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
#'     ymd("2020-12-06", "2020-01-01"),
#'     ymd("2020-11-11", NA)
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
  partialdate <- str_match(dtc, paste0(
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
  restricted <- restrict_imputed_dtc_dt(
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
#'
#' @family utils_impute
#'
#' @keywords internal
#'
#' @seealso [impute_dtc_dtm()], [impute_dtc_dt()]
restrict_imputed_dtc_dt <- function(dtc,
                                    imputed_dtc,
                                    min_dates,
                                    max_dates) {
  if (!(is.null(min_dates) || length(min_dates) == 0) ||
    !(is.null(max_dates) || length(max_dates) == 0)) {
    suppress_warning(
      { # nolint
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
      min_date_iso <- strftime(min_date, format = "%Y-%m-%d", tz = "UTC")
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
#' compute_dtf(dtc = "--06-01T00:00", dt = as.Date("2022-06-01"))
#' compute_dtf(dtc = "2022-06--T00:00", dt = as.Date("2022-06-01"))
#' compute_dtf(dtc = "2022---01T00:00", dt = as.Date("2022-06-01"))
#' compute_dtf(dtc = "2022----T00:00", dt = as.Date("2022-06-01"))
compute_dtf <- function(dtc, dt) {
  assert_character_vector(dtc)
  assert_date_vector(dt)

  is_na <- is.na(dt)
  n_chr <- nchar(dtc)
  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  # Find date portion
  date_portion <- ifelse(grepl("T", dtc),
    gsub("T", "", substr(dtc, 1, str_locate(dtc, "T")[, 1])),
    substr(dtc, 1, 10)
  )
  n_chr_date_portion <- nchar(date_portion)

  # Location of the first instance of the double hyphen to determine if its month/day imputation
  location_of_double_hyphen <- str_locate(date_portion, "--")[, 1]

  case_when(
    (!is_na & n_chr >= 10 & n_chr_date_portion == 10 & valid_dtc) | is_na | !valid_dtc ~ NA_character_, # nolint
    n_chr_date_portion < 4 | is.na(dtc) ~ "Y",
    n_chr_date_portion < 10 & location_of_double_hyphen == 1 ~ "Y", # dates like "--07-07"
    n_chr_date_portion == 4 ~ "M",
    n_chr_date_portion < 10 & location_of_double_hyphen == 5 ~ "M", # dates like "2019---07"
    n_chr_date_portion == 7 ~ "D",
    n_chr_date_portion < 10 & location_of_double_hyphen == 8 ~ "D", # dates like "2019-07--"
  )
}
