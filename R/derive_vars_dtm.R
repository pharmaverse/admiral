#' Derive/Impute a Datetime from a Character Date
#'
#' Derive a datetime object (`*DTM`) from a character date (`--DTC`).
#' The date and time can be imputed (see `date_imputation`/`time_imputation` arguments)
#' and the date/time imputation flag (`*DTF`, `*TMF`) can be added.
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
#'   the specified prefix, for the date imputation flag (`*DTF`), and for the time
#'   imputation flag (`*TMF`), i.e., for `new_vars_prefix = "AST"` the variables
#'   `ASTDT`, `ASTDTF`, and `ASTTMF` are created.
#'
#' @permitted [char_scalar]
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
#'  Please note that CDISC requirements dictate the need for a date/time imputation
#'  flag if any imputation is performed, so `flag_imputation = "none"` should
#'  only be used if the imputed variable is not part of the final ADaM dataset.
#'
#' @permitted [date_time_flag_imp]
#'
#' @inheritParams impute_dtc_dtm
#' @inheritParams compute_tmf
#'
#' @details
#' The presence of a `*DTF` variable is checked and the variable is not derived
#' if it already exists in the input dataset. However, if `*TMF` already exists
#' in the input dataset, a warning is issued and `*TMF` will be overwritten.
#'
#' @return  The input dataset with the datetime `*DTM` (and the date/time imputation
#' flag `*DTF`, `*TMF`) added.
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
#' @caption Derive a datetime variable imputing time
#' @info In this example, we derive `ASTDTM` from `MHSTDTC`. Note that by default the function
#' imputes missing time components to `00` but doesn't impute missing date components
#' and automatically produces the time imputation flag (`ASTTMF`).
#' @code
#' library(tibble)
#' library(lubridate)
#'
#' mhdt <- tribble(
#'   ~MHSTDTC,
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
#'   dtc = MHSTDTC
#' )
#'
#' @caption Impute to the latest (`date_imputation = "last"`)
#' @info In this example, we set `date_imputation = "last"` to get the last month/day
#' for partial dates. We also set `time_imputation = "last"`. The function will use
#' all or part of `23:59:59` for time imputation. Note that `highest_imputation` must
#' be at least `"D"` to perform date imputation. Here we use `highest_imputation = "M"`
#' to request imputation of month and day (and time). Also note that
#' two flag variables are created. By default `ASTTMF` is set to `NA`
#' if only seconds are imputed. Set `ignore_seconds_flag = FALSE`
#' to have the `"S"` flag for `ASTTMF`.
#'
#' @code
#' derive_vars_dtm(
#'  mhdt,
#'  new_vars_prefix = "AST",
#'  dtc = MHSTDTC,
#'  date_imputation = "last",
#'  time_imputation = "last",
#'  highest_imputation = "M"
#' )
#'
#'
#' @caption Suppress imputation flags (`flag_imputation = "none"`)
#' @info In this example, we derive `ASTDTM` but suppress the `ASTTMF`. Note that
#' function appends missing `"hh:mm:ss"` to `ASTDTM`. The `flag_imputation = "none"`
#' call ensures no date/time imputation flag is created. In practice, as per CDISC
#' requirements this option can only be selected if the imputed variable is not part
#' of the final ADaM dataset.
#' @code
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   flag_imputation = "none"
#' )
#'
#' @caption Avoid imputation after specified datetimes (`max_dates`)
#' @info In this example, we derive `AENDTM` where AE end date is imputed to the last date.
#' To ensure that the imputed date is not after the death or data cut off date we can
#' set `max_dates = exprs(DTHDT, DCUTDT)`. Note two flag variables: `ASTDTF` and `ASTTMF`
#' are created. Setting `highest_imputation = "Y"` will allow for the missing `AEENDTC`
#' record to be imputed from `max_dates = exprs(DTHDT, DCUTDT)`.
#' @code
#' adae <- tribble(
#'    ~AEENDTC,             ~DTHDT,           ~DCUTDT,
#'    "2020-12", ymd("2020-12-26"), ymd("2020-12-24"),
#'    "2020-11", ymd("2020-12-06"), ymd("2020-12-24"),
#'           "", ymd("2020-12-06"), ymd("2020-12-24"),
#' "2020-12-20", ymd("2020-12-06"), ymd("2020-12-24")
#' )
#'
#' derive_vars_dtm(
#'   adae,
#'   dtc = AEENDTC,
#'   new_vars_prefix = "AEN",
#'   highest_imputation = "Y",
#'   date_imputation = "last",
#'   time_imputation = "last",
#'   max_dates = exprs(DTHDT, DCUTDT)
#' )
#'
#' @caption Include `"S"` for time imputation flag (`ignore_seconds_flag`)
#' @info In this example, we set `ignore_seconds_flag = FALSE` to include `S` for
#' seconds in the `ASTTMF` variable. The default value of `ignore_seconds_flag`
#' is `TRUE` so the `"S"` is not normally displayed. The ADaM IG states that given
#' SDTM (`--DTC`) variable, if only hours and minutes are ever collected, and
#' seconds are imputed in (`*DTM`) as `00`, then it is not necessary to
#' set (`*TMF`) to `"S"`.
#' @code
#'
#' mhdt <- tribble(
#' ~MHSTDTC,
#' "2019-07-18T15:25",
#' "2019-07-18",
#' "2019-02",
#' "2019",
#' "2019---07",
#' ""
#' )
#'
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   ignore_seconds_flag = FALSE
#' )
#'
#' @caption Preserve lower components if higher ones were imputed (`preserve`)
#' @info In this example, we impute dates as the middle month/day with `date_imputation = "mid"`
#' and impute time as last (`23:59:59`) with `time_imputation = "last"`.
#' We use the `preserve` argument to "preserve" partial dates.  For example,
#' `"2019---18T15:-:05"`, will be displayed as `"2019-06-18 15:59:05"` by setting
#' `preserve = TRUE`.
#' @code
#' mhdt <- tribble(
#' ~MHSTDTC,
#' "2019-07-18T15:25",
#' "2019---18T15:-:05",
#' "2019-07-18",
#' "2019-02",
#' "2019",
#' "2019---07",
#' ""
#' )
#'
#' derive_vars_dtm(
#'   mhdt,
#'   new_vars_prefix = "AST",
#'   dtc = MHSTDTC,
#'   highest_imputation = "M",
#'   date_imputation = "mid",
#'   time_imputation = "last",
#'   preserve = TRUE,
#'   ignore_seconds_flag = FALSE
#' )
#' @caption Further examples
#' @info Further example usages of this function can be found in the
#'   `vignette("imputation")`.
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
                            ignore_seconds_flag = TRUE) {
  # check and quote arguments
  dtc <- assert_symbol(enexpr(dtc))
  assert_data_frame(dataset, required_vars = exprs(!!dtc))

  assert_character_scalar(new_vars_prefix)

  flag_imputation <- assert_character_scalar(
    flag_imputation,
    values = c("auto", "both", "date", "time", "none"),
    case_sensitive = FALSE
  )

  # the `assert_dt_dtm_inputs` function is stored in `derive_vars_dt_dtm_utils.R`
  assert_highest_imputation(
    highest_imputation = highest_imputation,
    highest_imputation_values = c("Y", "M", "D", "h", "m", "s", "n"),
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates
  )

  dtm <- paste0(new_vars_prefix, "DTM")

  # Issue a warning if *DTM already exists
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
    # add *DTF if not there already
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
    # add *TMF variable
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
#' Convert a date character vector (usually `--DTC`) into a Date vector (usually `*DTM`).
#'
#' @param dtc The `--DTC` date to convert.
#'
#' @permitted [date_chr_vector]
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

#' Impute Partial Date(-time) Portion of a `--DTC` Variable
#'
#' Imputation partial date/time portion of a `--DTC` variable. based on user
#' input.
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
#'   If `"n"` is specified, no imputation is performed, i.e., if any component is
#'   missing, `NA_character_` is returned.
#'
#'   If `"Y"` is specified, `date_imputation` should be `"first"` or `"last"`
#'   and `min_dates` or `max_dates` should be specified respectively. Otherwise,
#'   `NA_character_` is returned if the year component is missing.
#'
#' @permitted [date_time_high_imp]
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
#' @permitted [time_imp]
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
#'    ymd_hm("2020-12-06T12:12"),
#'    ymd_hm("2020-11-11T11:11")
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
#' @permitted [date_list]
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
#'
#' @permitted [date_list]
#'
#' @param preserve Preserve lower level date/time part when higher order part
#' is missing, e.g. preserve day if month is missing or
#' preserve minute when hour is missing.
#'
#' For example `"2019---07"` would return `"2019-06-07` if `preserve = TRUE`
#' (and `date_imputation = "mid"`).
#'
#' @permitted [boolean]
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
#'     ymd_hm("2020-12-06T12:12"),
#'     ymd_hm("2020-11-11T11:11")
#'   ),
#'   highest_imputation = "M"
#' )
#'
#' # Impute completely missing dates (only possible if min_dates or max_dates is specified)
#' impute_dtc_dtm(
#'   c("2020-12", NA_character_),
#'   min_dates = list(
#'     ymd_hm("2020-12-06T12:12", "2020-01-01T01:01"),
#'     ymd_hm("2020-11-11T11:11", NA)
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

  assert_highest_imputation(
    highest_imputation = highest_imputation,
    highest_imputation_values = imputation_levels,
    date_imputation = date_imputation,
    min_dates = min_dates,
    max_dates = max_dates
  )


  highest_imputation <- dtm_level(highest_imputation)

  # the `assert_date_imputation` function is stored in `derive_vars_dt_dtm_utils.R`
  date_imputation <- assert_date_imputation(
    highest_imputation = highest_imputation,
    date_imputation = date_imputation
  )

  # the `assert_time_imputation` function is stored in `derive_vars_dt_dtm_utils.R`
  time_imputation <- assert_time_imputation(
    highest_imputation = highest_imputation,
    time_imputation = time_imputation
  )
  assert_logical_scalar(preserve)

  if (length(dtc) == 0) {
    return(vector("character"))
  }

  # Parse character date ----
  partial <- get_partialdatetime(dtc, create_datetime = TRUE)
  components <- names(partial)

  # Handle preserve argument ----
  if (!preserve) {
    partial <- propagate_na_values(partial)
  }

  # Determine target components ----
  target <- get_imputation_targets(partial, date_imputation, time_imputation)

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

  imputed <- impute_date_time(partial, target)
  imputed_dtc <- format_imputed_dtc(imputed)

  if (date_imputation == "last") {
    imputed_dtc <- adjust_last_day_imputation(imputed_dtc, partial)
  }

  # Handle min_dates and max_dates argument ----
  restricted <- restrict_imputed_dtc_dtm(
    dtc,
    imputed_dtc = imputed_dtc,
    min_dates = min_dates,
    max_dates = max_dates
  )

  restricted
}

#' Restrict Imputed `--DTC` date to Minimum/Maximum Dates
#'
#' @param imputed_dtc The imputed `--DTC` date
#'
#' @inheritParams impute_dtc_dtm
#'
#' @returns
#'   - The last of the minimum dates (`min_dates`) which are in the range of the
#'   partial `--DTC` date (`dtc`)
#'   - The first of the maximum dates (`max_dates`) which are in the range of the
#'   partial `--DTC` date (`dtc`)
#'   - `imputed_dtc` if the partial `--DTC` date (`dtc`) is not in range of any of
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
  any_mindate <- !(is.null(min_dates) || length(min_dates) == 0)
  any_maxdate <- !(is.null(max_dates) || length(max_dates) == 0)
  if (any_mindate || any_maxdate) {
    dtc_range <-
      get_dt_dtm_range(
        dtc,
        create_datetime = TRUE
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
      min_date_iso <- strftime(min_date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
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
#' Derive the time imputation flag (`*TMF`) comparing a date character vector
#' (`--DTC`) with a Datetime vector (`*DTM`).
#'
#' @param dtc The date character vector (`--DTC`).
#'
#'   A character date is expected in a format like `yyyy-mm-ddThh:mm:ss` (partial or complete).
#'
#' @param dtm The Date vector to compare (`*DTM`).
#'
#'   A datetime object is expected.
#'
#' @param ignore_seconds_flag  ADaM IG states that given SDTM (`--DTC`) variable,
#' if only hours and minutes are ever collected, and seconds are imputed in
#' (`*DTM`) as 00, then it is not necessary to set (`*TMF`) to `"S"`.
#'
#' By default it is assumed that no seconds are collected and `*TMF` shouldn't be set to `"S"`.
#' A user can set this to `FALSE` if seconds are collected.
#'
#' The default value of `ignore_seconds_flag` is set to `TRUE` in
#' admiral 1.4.0 and later.
#'
#' @permitted [boolean]
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return The time imputation flag (`*TMF`) (character value of `"H"`, `"M"` , `"S"` or `NA`)
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
#' compute_tmf(dtc = "2019-07-18T15:25", dtm = ymd_hm("2019-07-18T15:25"))
#' compute_tmf(dtc = "2019-07-18T15", dtm = ymd_hm("2019-07-18T15:25"))
#' compute_tmf(dtc = "2019-07-18", dtm = ymd("2019-07-18"))
#' compute_tmf(dtc = "2022-05--T00:00", dtm = ymd_hm("2022-05-15T23:59"))
#' compute_tmf(dtc = "2022-05--T23:00", dtm = ymd_hm("2022-05-15T23:59"))
#' compute_tmf(
#'   dtc = "2022-05--T23:59:00",
#'   dtm = ymd_hms("2022-05-15T23:59:59"),
#'   ignore_seconds_flag = FALSE
#' )
#'
compute_tmf <- function(dtc,
                        dtm,
                        ignore_seconds_flag = TRUE) {
  assert_date_vector(dtm)
  assert_character_vector(dtc)
  assert_logical_scalar(ignore_seconds_flag)

  valid_dtc <- is_valid_dtc(dtc)
  warn_if_invalid_dtc(dtc, valid_dtc)

  partial <- get_partialdatetime(dtc, create_datetime = TRUE)
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
