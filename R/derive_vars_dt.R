#' Derive/Impute a Date from a Character Date
#'
#' Derive a date (`*DT`) from a character date (`--DTC`).
#' The date can be imputed (see `date_imputation` argument)
#' and the date imputation flag (`*DTF`) can be added.
#'
#' In `{admiral}` we don't allow users to pick any single part of the date/time to
#' impute, we only enable to impute up to a highest level, i.e. you couldn't
#' choose to say impute months, but not days.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("dtc"))`
#'
#' @permitted [dataset]
#'
#' @param new_vars_prefix Prefix used for the output variable(s).
#'
#'   A character scalar is expected. For the date variable (`*DT`) is appended to
#'   the specified prefix and for the date imputation flag (`*DTF`), i.e., for
#'   `new_vars_prefix = "AST"` the variables `ASTDT` and `ASTDTF` are created.
#'
#' @permitted [char_scalar]
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
#'  Please note that CDISC requirements dictate the need for a date imputation
#'  flag if any imputation is performed, so `flag_imputation = "none"` should
#'  only be used if the imputed variable is not part of the final ADaM dataset.
#'
#' @permitted [date_flag_imp]
#'
#' @inheritParams impute_dtc_dt
#'
#' @return
#' The input dataset with the date `*DT` (and the date imputation flag `*DTF`
#' if requested) added.
#'
#' @details
#' The presence of a `*DTF` variable is checked and if it already exists in the input dataset,
#' a warning is issued and `*DTF` will be overwritten.
#'
#'
#' @seealso `vignette("imputation")`
#'
#' @family der_date_time
#'
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @examplesx
#'
#' @caption Derive a date variable without imputation
#'
#' @info In this example, we derive `ASTDT` from `MHSTDTC` with no imputation
#' done for partial dates.
#' @code
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
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC
#' )
#'
#' @caption Impute partial dates (`highest_imputation`)
#' @info Imputation is requested by the `highest_imputation` argument. Here
#' `highest_imputation = "M"` for month imputation is used, i.e. the highest
#' imputation done on a partial date is up to the month. By default, missing date
#' components are imputed to the first day/month/year. A date imputation flag variable, `ASTDTF`,
#' is automatically created. The flag variable indicates if imputation was done
#' on the date.
#'
#' @code
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "first"
#' )
#'
#' @caption Impute to the last day/month (`date_imputation = "last"`)
#' @info In this example, we derive `ADT` impute partial dates to last day/month, i.e.
#' `date_imputation = "last"`.
#'
#' @code
#' qsdt <- tribble(
#'   ~QSDTC,
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#'
#' derive_vars_dt(
#'   qsdt,
#'   new_vars_prefix = "A",
#'   dtc = QSDTC,
#'   highest_imputation = "M",
#'   date_imputation = "last"
#' )
#'
#' @caption Impute to the middle (`date_imputaton = "mid"`) and suppress
#' imputation flag (`flag_imputation = "none"`)
#' @info In this example, we will derive `TRTSDT` with date imputation flag
#' (`*DTF`) suppressed. Since `date_imputation = "mid"`, partial date imputation
#'  will be set to June 30th for missing month and 15th for missing day only.
#'  The `flag_imputation = "none"` call ensures no date imputation flag is
#'  created. In practice, as per CDISC requirements this option can only be
#'  selected if the imputed variable is not part of the final ADaM dataset.
#'
#' @code
#'
#' exdt <- tribble(
#'   ~EXSTDTC,
#'   "2019-07-18T15:25:40",
#'   "2019-07-18T15:25",
#'   "2019-07-18",
#'   "2019-02",
#'   "2019",
#'   "2019---07",
#'   ""
#' )
#' derive_vars_dt(
#'   exdt,
#'   new_vars_prefix = "TRTS",
#'   dtc = EXSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   flag_imputation = "none"
#' )
#'
#' @caption Impute to a specific date (`date_imputation = "04-06"`)
#' @info In this example, we derive `ASTDT` with specific date imputation, i.e.
#' `date_imputation = "04-06"`. Note that day portion, `"-06"`, is used in the
#' imputation of the record with `"2019-02"`.
#'
#' @code
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "04-06"
#' )
#' @caption Avoid imputation before a user-defined date (`min_dates`)
#' @info In this example, we derive `ASTDT` where `AESTDTC` is all partial dates in
#' need of imputation. Using `min_dates = exprs(TRTSDTM)`, we are telling the function
#' to not allow imputation dates to be before the treatment start date
#' via `min_dates` argument. Note that the second record does not get imputed
#' as it is before `TRTSDTM`.
#'
#' @code
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
#' @caption Preserve lower components if higher ones were imputed (`preserve`)
#' @info The `preserve` argument can be used to "preserve" information from the partial dates.
#' For example, `"2019---07"`, will be displayed as `"2019-06-07"` rather than
#' `"2019-06-30"` with `preserve = TRUE` and `date_imputation = "mid"` .
#'
#' @code
#' derive_vars_dt(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   preserve = TRUE
#' )
#' @caption Further examples
#' @info Further example usages of this function can be found in the
#'   `vignette("imputation")`.
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
  dtc <- assert_symbol(enexpr(dtc))
  assert_data_frame(dataset, required_vars = exprs(!!dtc))

  assert_character_scalar(new_vars_prefix)

  flag_imputation <- assert_character_scalar(
    flag_imputation,
    values = c("auto", "date", "none"),
    case_sensitive = FALSE
  )

  # the `assert_highest_imputation` function is stored in `derive_vars_dt_dtm_utils.R`
  assert_highest_imputation(
    highest_imputation = highest_imputation,
    highest_imputation_values = c("Y", "M", "D", "n"),
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates
  )

  # output varname
  dt <- paste0(new_vars_prefix, "DT")
  warn_if_vars_exist(dataset, dt)

  # derive *DT var
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
    # add *DTF if not there already
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
#' Convert a date character vector (usually `--DTC`) into a Date vector (usually `*DT`).
#'
#' @param dtc The --DTC date to convert.
#'
#' @permitted [date_chr_vector]
#'
#' @inheritParams impute_dtc_dt
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

#' Impute Partial Date Portion of a `--DTC` Variable
#'
#' Imputation partial date portion of a `--DTC` variable based on user input.
#'
#' @param dtc The `--DTC` date to impute
#'
#'   A character date is expected in a format like `yyyy-mm-dd` or
#'   `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
#'   valid "missing" value for any component.
#'
#' @permitted [date_chr]
#'
#' @param highest_imputation Highest imputation level
#'
#'   The `highest_imputation` argument controls which components of the `--DTC`
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
#'   If `"Y"` is specified, `date_imputation` must be `"first"` or `"last"`
#'   and `min_dates` or `max_dates` must be specified respectively. Otherwise,
#'   an error is thrown.
#'
#' @permitted [date_high_imp]
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   A character value is expected.
#'    - If  `highest_imputation` is `"M"`, month and day can be
#'      specified as `"mm-dd"`: e.g. `"06-15"` for the 15th of June
#'    - When  `highest_imputation` is `"M"` or  `"D"`, the following keywords are available:
#'      `"first"`, `"mid"`, `"last"` to impute to the first/mid/last
#'      day/month. If `"mid"` is specified, missing components are imputed as the
#'      middle of the possible range:
#'       - If both month and day are missing, they are imputed as `"06-30"`
#'        (middle of the year).
#'       - If only day is missing, it is imputed as `"15"` (middle of the month).
#'
#'   The year can not be specified; for imputing the year
#'   `"first"` or `"last"` together with `min_dates` or `max_dates` argument can
#'   be used (see examples).
#'
#' @permitted [date_imp]
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
#' @permitted [date_list]
#'
#' @param max_dates Maximum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not after
#' any of the specified dates, e.g., that the imputed date is not after the data
#' cut off date. Only dates which are in the range of possible dates are
#' considered. A date or date-time object is expected.
#'
#' @permitted [date_list]
#'
#' @param preserve Preserve day if month is missing and day is present
#'
#' For example `"2019---07"` would return `"2019-06-07` if `preserve = TRUE`
#' (and `date_imputation = "MID"`).
#'
#' @permitted [boolean]
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

  assert_highest_imputation(
    highest_imputation = highest_imputation,
    highest_imputation_values = imputation_levels,
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates
  )

  highest_imputation <- dt_level(highest_imputation)

  assert_logical_scalar(preserve)

  # the `assert_date_imputation` function is stored in `derive_vars_dt_dtm_utils.R`
  date_imputation <- assert_date_imputation(
    highest_imputation = highest_imputation,
    date_imputation = date_imputation
  )

  if (length(dtc) == 0) {
    return(character(0))
  }

  # Parse partials
  partial <- get_partialdatetime(dtc, create_datetime = FALSE)
  components <- names(partial)

  # Handle preserve argument ----
  if (!preserve) {
    partial <- propagate_na_values(partial)
  }

  # Determine target components ----
  target <- get_imputation_targets(partial, date_imputation)

  for (c in components) {
    if (highest_imputation < dt_level(imputation_levels[[c]])) {
      target[[c]] <- "xx"
    }
  }

  # Impute ----
  imputed <- impute_date_time(partial, target)
  imputed_dtc <- format_imputed_dtc(imputed)


  if (date_imputation == "last") {
    imputed_dtc <- adjust_last_day_imputation(imputed_dtc, partial)
  }

  # Handle min_dates and max_dates argument ----
  restricted <- restrict_imputed_dtc_dt(
    dtc,
    imputed_dtc = imputed_dtc,
    min_dates = min_dates,
    max_dates = max_dates
  )

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
  any_mindate <- !(is.null(min_dates) || length(min_dates) == 0)
  any_maxdate <- !(is.null(max_dates) || length(max_dates) == 0)
  if (any_mindate || any_maxdate) {
    # determine range of possible dates
    dtc_range <-
      get_dt_dtm_range(
        dtc,
        create_datetime = FALSE
      )
    min_dtc <- dtc_range[["lower"]]
    max_dtc <- dtc_range[["upper"]]
  }
  if (any_mindate) {
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
  if (any_maxdate) {
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
#' Derive the date imputation flag (`*DTF`) comparing a date character vector
#' (`--DTC`) with a Date vector (`*DT`).
#'
#' @param dtc The date character vector (`--DTC`).
#'
#'   A character date is expected in a format like `yyyy-mm-ddThh:mm:ss` (partial or complete).
#'
#' @param dt The  Date vector to compare.
#'
#'   A date object is expected.
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The date imputation flag (`*DTF`) (character value of `"D"`, `"M"` , `"Y"` or `NA`)
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
