#' Convert `XXTPT` Strings to Hours
#'
#' Converts CDISC timepoint strings (e.g., `PCTPT`, `VSTPT`, `EGTPT`,
#' `ISTPT`, `LBTPT`) into numeric hours for analysis. The function handles
#' common dose-centric formats including pre-dose, post-dose (hours/minutes),
#' days, time ranges, and infusion-related time markers.
#'
#' @param xxtpt A character vector of timepoint descriptions from SDTM
#'   `--TPT` variables (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`).
#'   Can contain `NA` values.
#'
#' @permitted character vector
#'
#' @param infusion_duration Numeric value specifying the duration of infusion
#'   in hours. Used to convert "EOI" (End of Infusion) patterns and patterns
#'   describing time after end of infusion. Must be a positive number.
#'   Default is 1 hour.
#'
#' @permitted numeric scalar (positive)
#'
#' @param range_method Method for converting time ranges to single values.
#'   Options are "midpoint" (default), "start", or "end". For example,
#'   "0-6h" with midpoint returns 3, with start returns 0, with end returns 6.
#'
#' @permitted character scalar ("midpoint", "start", or "end")
#'
#' @details
#' The function recognizes the following patterns (all case-insensitive):
#'
#' **Special Cases:**
#' * `"Screening"` → -1
#' * `"Pre-dose"`, `"Predose"`, `"Pre-infusion"`, `"Pre-inf"`, `"Infusion"`,
#'   `"0H"` → 0
#' * `"EOI"`, `"End of Infusion"`, `"After End of Infusion"`,
#'   `"After End of Treatment"` → `infusion_duration` (default: 1)
#' * `"Morning"`, `"Evening"` → `NA_real_`
#' * Unrecognized values → `NA_real_`
#'
#' **Time Ranges:**
#' Time ranges are converted based on the `range_method` parameter:
#' * `"0-6h Post-dose"` with `range_method = "midpoint"` (default) → 3
#' * `"0-6h Post-dose"` with `range_method = "start"` → 0
#' * `"0-6h Post-dose"` with `range_method = "end"` → 6
#' * `"0-4H PRIOR START OF INFUSION"` with midpoint → -2 (negative for prior)
#' * `"8-16H POST START OF INFUSION"` with midpoint → 12
#'
#' **Time-based Conversions:**
#' * **Days**: `"Day 1"` → 24, `"30 DAYS AFTER LAST"` → 720
#' * **Hours + Minutes**: `"1H30M"` → 1.5
#' * **Hours**: `"2 hours"` → 2, `"1 HOUR POST"` → 1
#' * **Minutes**: `"30M"` → 0.5, `"30 MIN POST"` → 0.5
#' * **Predose**: `"5 MIN PREDOSE"` → -0.0833
#' * **Post EOI/EOT**: `"1 HOUR POST EOI"` → infusion_duration + 1,
#'   `"24 HR POST INF"` → infusion_duration + 24
#' * **After end**: `"30MIN AFTER END OF INFUSION"` → infusion_duration + 0.5
#' * **Start of infusion**: `"8H PRIOR START OF INFUSION"` → -8
#' * **Pre EOI**: `"10MIN PRE EOI"` → -0.1667
#'
#' **Supported Unit Formats:**
#' * Hours: H, h, HR, hr, HOUR, hour (with optional plurals)
#' * Minutes: M, m, MIN, min, MINUTE, minute (with optional plurals)
#' * Days: D, d, DAY, day (with optional plurals)
#' * Flexible whitespace and optional "Post-dose", "POST", "After last"
#'   suffixes
#'
#' @return A numeric vector of timepoints in hours. Returns `NA_real_` for:
#' * Input `NA` values
#' * Unrecognized timepoint formats
#' * Non-time descriptors (e.g., "Morning", "Evening")
#'
#' Returns `numeric(0)` for empty input.
#'
#' @keywords com_date_time
#' @family com_date_time
#'
#' @export
#'
#' @examples
#' # Basic dose-centric patterns
#' convert_xxtpt_to_hours(c(
#'   "Screening",
#'   "Pre-dose",
#'   "Pre-inf",
#'   "30M",
#'   "1H",
#'   "2H POSTDOSE",
#'   "Day 1"
#' ))
#'
#' # Predose patterns (negative times)
#' convert_xxtpt_to_hours(c("5 MIN PREDOSE", "1 HOUR PREDOSE"))
#'
#' # Infusion-related patterns
#' convert_xxtpt_to_hours(c(
#'   "1 HOUR POST EOI",
#'   "24 HR POST INF",
#'   "30MIN AFTER END OF INFUSION",
#'   "8H PRIOR START OF INFUSION",
#'   "10MIN PRE EOI"
#' ))
#'
#' # Time ranges - default midpoint
#' convert_xxtpt_to_hours(c(
#'   "0-6h Post-dose",
#'   "0-4H PRIOR START OF INFUSION",
#'   "8-16H POST START OF INFUSION"
#' ))
#'
#' # Time ranges - specify method
#' convert_xxtpt_to_hours("0-6h Post-dose", range_method = "end")
#' convert_xxtpt_to_hours("0-6h Post-dose", range_method = "start")
#'
#' # Custom infusion duration (2 hours)
#' convert_xxtpt_to_hours(
#'   c("EOI", "End of Infusion", "1 HOUR POST EOI"),
#'   infusion_duration = 2
#' )
convert_xxtpt_to_hours <- function(xxtpt,
                                   infusion_duration = 1,
                                   range_method = "midpoint") {
  # Validate inputs
  assert_character_vector(xxtpt)
  assert_numeric_vector(infusion_duration, length = 1)
  assert_character_vector(range_method, values = c("start", "end", "midpoint"))

  if (infusion_duration <= 0) {
    cli_abort(
      "{.arg infusion_duration} must be positive, but is {infusion_duration}."
    )
  }

  if (length(xxtpt) == 0) {
    return(numeric(0))
  }

  xxtpt <- trimws(xxtpt)
  result <- rep(NA_real_, length(xxtpt))
  na_idx <- is.na(xxtpt)

  # Process patterns in order of specificity
  result <- convert_special_cases(xxtpt, result, na_idx, infusion_duration)
  result <- convert_time_units(xxtpt, result, na_idx)
  result <- convert_ranges(xxtpt, result, na_idx, range_method)
  result <- convert_infusion_patterns(xxtpt, result, na_idx, infusion_duration)
  result <- convert_simple_units(xxtpt, result, na_idx)

  result
}

# Helper Functions ----

#' Convert Special Case Patterns
#'
#' Converts special case timepoint patterns to numeric hours including screening,
#' pre-dose, pre-infusion, and end of infusion markers.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param infusion_duration Duration of infusion in hours (positive numeric)
#'
#' @details
#' Recognizes and converts the following patterns:
#' * "Screening" → -1
#' * "Pre-dose", "Predose", "Pre-infusion", "Pre-inf", "Infusion", "0H" → 0
#' * "EOI", "End of Infusion", "After End of Infusion/Treatment" → infusion_duration
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with special case patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_special_cases <- function(xxtpt, result, na_idx, infusion_duration) {
  # Screening
  screening_pattern <- regex("^screening$", ignore_case = TRUE)
  screening_idx <- str_detect(xxtpt, screening_pattern) & !na_idx
  result[screening_idx] <- -1

  # Pre-dose, Predose, Pre-infusion, Pre-inf, Infusion, 0H -> 0
  zero_pattern <- regex(
    "^(pre-?dose|pre-?inf(?:usion)?|infusion|0\\s*h(?:r|our)?s?)$",
    ignore_case = TRUE
  )
  zero_idx <- str_detect(xxtpt, zero_pattern) & is.na(result) & !na_idx
  result[zero_idx] <- 0

  # EOI, End of Infusion, After End of Infusion/Treatment -> infusion_duration
  eoi_pattern <- regex(
    paste0(
      "^(eoi|end\\s+of\\s+infusion|",
      "after\\s+end\\s+of\\s+(?:infusion|treatment))$"
    ),
    ignore_case = TRUE
  )
  eoi_idx <- str_detect(xxtpt, eoi_pattern) & is.na(result) & !na_idx
  result[eoi_idx] <- infusion_duration

  result
}

#' Convert Time Unit Patterns
#'
#' Converts days and hours+minutes combination patterns to numeric hours.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes and converts the following patterns:
#' * Days: "Day 1" → 24, "2D" → 48, "30 DAYS AFTER LAST" → 720
#' * Hours+minutes: "1H30M" → 1.5, "2HR15MIN" → 2.25
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with time unit patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_time_units <- function(xxtpt, result, na_idx) {
  # Days (convert to hours by multiplying by 24)
  days_pattern <- regex(
    paste0(
      "^(?:day\\s+)?(?<days>\\d+(?:\\.\\d+)?)\\s*",
      "(?:d|day|days)?",
      "(?:\\s+(?:after\\s+last|post(?:\\s*-?\\s*dose)?))?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  days_matches <- str_match(xxtpt, days_pattern)
  days_idx <- !is.na(days_matches[, 1]) & is.na(result) & !na_idx
  if (any(days_idx)) {
    result[days_idx] <- as.numeric(days_matches[days_idx, "days"]) * 24
  }

  # Hours+minutes combinations
  hm_pattern <- regex(
    paste0(
      "^(?<hours>\\d+(?:\\.\\d+)?)\\s*h(?:r|our)?s?\\s*",
      "(?<minutes>\\d+(?:\\.\\d+)?)\\s*m(?:in|inute)?s?",
      "(?:\\s+post(?:\\s*-?\\s*dose)?)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hm_matches <- str_match(xxtpt, hm_pattern)
  hm_idx <- !is.na(hm_matches[, 1]) & is.na(result) & !na_idx
  if (any(hm_idx)) {
    hours <- as.numeric(hm_matches[hm_idx, "hours"])
    minutes <- as.numeric(hm_matches[hm_idx, "minutes"])
    result[hm_idx] <- hours + minutes / 60
  }

  result
}

#' Convert Range Patterns
#'
#' Converts time range patterns to numeric hours using specified range method
#' (start, midpoint, or end).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param range_method Method for selecting value from range: "start", "midpoint",
#'   or "end"
#'
#' @details
#' Recognizes and converts the following range patterns:
#' * Ranges with direction: "0-4H PRIOR START OF INFUSION" → -2 (with midpoint)
#' * Simple ranges: "0-6h Post-dose" → 3 (with midpoint)
#'
#' Direction affects sign: "PRIOR" results in negative values, "POST" in positive.
#' Processes ranges with direction before simple ranges to ensure correct matching.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with range patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_ranges <- function(xxtpt, result, na_idx, range_method) {
  # Ranges with direction (PRIOR/POST START/END)
  range_dir_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?\\s+(?<direction>prior|post)\\s+(?:start|end)"
    ),
    ignore_case = TRUE
  )
  range_dir_matches <- str_match(xxtpt, range_dir_pattern)
  range_dir_idx <- !is.na(range_dir_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_dir_idx)) {
    start_val <- as.numeric(range_dir_matches[range_dir_idx, "start"])
    end_val <- as.numeric(range_dir_matches[range_dir_idx, "end"])
    direction <- tolower(range_dir_matches[range_dir_idx, "direction"])

    range_val <- calculate_range_value(start_val, end_val, range_method)
    result[range_dir_idx] <- if_else(
      direction == "prior",
      -range_val,
      range_val
    )
  }

  # Simple time ranges
  range_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?(?:\\s+post(?:\\s*-?\\s*dose)?)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  range_matches <- str_match(xxtpt, range_pattern)
  range_idx <- !is.na(range_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_idx)) {
    start_val <- as.numeric(range_matches[range_idx, "start"])
    end_val <- as.numeric(range_matches[range_idx, "end"])
    result[range_idx] <- calculate_range_value(start_val, end_val, range_method)
  }

  result
}

#' Calculate Range Value Based on Method
#'
#' Calculates a single numeric value from a range (start and end) using the
#' specified method.
#'
#' @param start_val Numeric vector of start values
#' @param end_val Numeric vector of end values (same length as start_val)
#' @param range_method Method for calculation: "start", "midpoint", or "end"
#'
#' @details
#' * "start": Returns start_val
#' * "end": Returns end_val
#' * "midpoint": Returns (start_val + end_val) / 2
#'
#' @return Numeric vector of calculated values
#'
#' @keywords internal
#' @noRd
calculate_range_value <- function(start_val, end_val, range_method) {
  switch(range_method,
    start = start_val,
    end = end_val,
    midpoint = (start_val + end_val) / 2
  )
}

#' Convert Infusion-Related Patterns
#'
#' Orchestrates conversion of all infusion-related timepoint patterns by calling
#' specialized helper functions for each pattern type.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param infusion_duration Duration of infusion in hours (positive numeric)
#'
#' @details
#' Calls helper functions in sequence to convert:
#' * Predose patterns (negative times)
#' * Post EOI/EOT patterns
#' * After end of infusion/treatment patterns
#' * HR POST INF patterns
#' * Start of infusion patterns
#' * MIN AFTER START patterns
#' * MIN PRE EOI patterns
#'
#' @return Updated numeric vector with all infusion-related patterns converted
#'
#' @keywords internal
#' @noRd
convert_infusion_patterns <- function(xxtpt, result, na_idx, infusion_duration) {
  result <- convert_predose_patterns(xxtpt, result, na_idx)
  result <- convert_post_eoi_patterns(xxtpt, result, na_idx, infusion_duration)
  result <- convert_after_end_patterns(xxtpt, result, na_idx, infusion_duration)
  result <- convert_hr_post_patterns(xxtpt, result, na_idx, infusion_duration)
  result <- convert_start_inf_patterns(xxtpt, result, na_idx)
  result <- convert_min_after_start(xxtpt, result, na_idx)
  result <- convert_min_pre_eoi(xxtpt, result, na_idx)
  result
}

#' Convert Predose Patterns
#'
#' Converts predose timepoint patterns to negative numeric hours (time before dose).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "5 MIN PREDOSE" → -0.0833 (negative 5 minutes)
#' * "1 HOUR PREDOSE" → -1 (negative 1 hour)
#'
#' Returns negative values to indicate time before dose.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with predose patterns converted to negative hours
#'
#' @keywords internal
#' @noRd
convert_predose_patterns <- function(xxtpt, result, na_idx) {
  predose_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+predose$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  predose_matches <- str_match(xxtpt, predose_pattern)
  predose_idx <- !is.na(predose_matches[, 1]) & is.na(result) & !na_idx
  if (any(predose_idx)) {
    time_value <- as.numeric(predose_matches[predose_idx, "value"])
    unit <- tolower(predose_matches[predose_idx, "unit"])
    is_minutes <- substr(unit, 1, 1) == "m"
    result[predose_idx] <- if_else(
      is_minutes,
      -time_value / 60,
      -time_value
    )
  }
  result
}

#' Convert Post EOI/EOT Patterns
#'
#' Converts "post end of infusion/treatment" patterns to numeric hours, adding
#' the infusion duration to account for time after infusion completion.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param infusion_duration Duration of infusion in hours (positive numeric)
#'
#' @details
#' Recognizes patterns like:
#' * "1 HOUR POST EOI" → infusion_duration + 1
#' * "30 MIN POST EOT" → infusion_duration + 0.5
#'
#' Adds infusion_duration because these patterns describe time AFTER the end of
#' infusion, so total time from dose start includes the infusion duration plus
#' the additional time.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with post EOI/EOT patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_post_eoi_patterns <- function(xxtpt, result, na_idx, infusion_duration) {
  post_eoi_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+post\\s+(?:eoi|eot)$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  post_eoi_matches <- str_match(xxtpt, post_eoi_pattern)
  post_eoi_idx <- !is.na(post_eoi_matches[, 1]) & is.na(result) & !na_idx
  if (any(post_eoi_idx)) {
    time_value <- as.numeric(post_eoi_matches[post_eoi_idx, "value"])
    unit <- tolower(post_eoi_matches[post_eoi_idx, "unit"])
    is_minutes <- substr(unit, 1, 1) == "m"
    result[post_eoi_idx] <- if_else(
      is_minutes,
      infusion_duration + time_value / 60,
      infusion_duration + time_value
    )
  }
  result
}

#' Convert After End of Infusion/Treatment Patterns
#'
#' Converts "after end of infusion/treatment" patterns (compact format) to numeric
#' hours, adding the infusion duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param infusion_duration Duration of infusion in hours (positive numeric)
#'
#' @details
#' Recognizes patterns like:
#' * "30MIN AFTER END OF INFUSION" → infusion_duration + 0.5
#' * "1 HOUR AFTER END OF TREATMENT" → infusion_duration + 1
#'
#' Adds infusion_duration because these patterns describe time AFTER the end of
#' infusion.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with after end patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_after_end_patterns <- function(xxtpt, result, na_idx, infusion_duration) {
  after_end_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "after\\s+end\\s+of\\s+(?:infusion|treatment)$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  after_end_matches <- str_match(xxtpt, after_end_pattern)
  after_end_idx <- !is.na(after_end_matches[, 1]) & is.na(result) & !na_idx
  if (any(after_end_idx)) {
    time_value <- as.numeric(after_end_matches[after_end_idx, "value"])
    unit <- tolower(after_end_matches[after_end_idx, "unit"])
    is_minutes <- substr(unit, 1, 1) == "m"
    result[after_end_idx] <- if_else(
      is_minutes,
      infusion_duration + time_value / 60,
      infusion_duration + time_value
    )
  }
  result
}

#' Convert HR POST INF/EOI/EOT Patterns
#'
#' Converts "HR POST INF/EOI/EOT" patterns to numeric hours, adding the infusion
#' duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param infusion_duration Duration of infusion in hours (positive numeric)
#'
#' @details
#' Recognizes patterns like:
#' * "24 HR POST INF" → infusion_duration + 24
#' * "24 HR POST EOI" → infusion_duration + 24
#'
#' Adds infusion_duration because these patterns describe time AFTER the end of
#' infusion.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with HR POST patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_hr_post_patterns <- function(xxtpt, result, na_idx, infusion_duration) {
  hr_post_pattern <- regex(
    "^(?<value>\\d+(?:\\.\\d+)?)\\s*hr\\s+post\\s+(?:inf|eoi|eot)$",
    ignore_case = TRUE,
    comments = TRUE
  )
  hr_post_matches <- str_match(xxtpt, hr_post_pattern)
  hr_post_idx <- !is.na(hr_post_matches[, 1]) & is.na(result) & !na_idx
  if (any(hr_post_idx)) {
    result[hr_post_idx] <- infusion_duration +
      as.numeric(hr_post_matches[hr_post_idx, "value"])
  }
  result
}

#' Convert Start of Infusion Patterns
#'
#' Converts patterns relative to the start of infusion to numeric hours, with
#' sign based on direction (prior = negative, post = positive).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "8H PRIOR START OF INFUSION" → -8
#' * "8H POST START OF INFUSION" → 8
#'
#' Direction affects sign: "PRIOR" results in negative values, "POST" in positive.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with start of infusion patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_start_inf_patterns <- function(xxtpt, result, na_idx) {
  h_start_inf_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)h\\s+",
      "(?<direction>prior|post)\\s+start\\s+of\\s+infusion$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  h_start_inf_matches <- str_match(xxtpt, h_start_inf_pattern)
  h_start_inf_idx <- !is.na(h_start_inf_matches[, 1]) & is.na(result) & !na_idx
  if (any(h_start_inf_idx)) {
    time_value <- as.numeric(h_start_inf_matches[h_start_inf_idx, "value"])
    direction <- tolower(h_start_inf_matches[h_start_inf_idx, "direction"])
    result[h_start_inf_idx] <- if_else(
      direction == "prior",
      -time_value,
      time_value
    )
  }
  result
}

#' Convert MIN AFTER START INF Patterns
#'
#' Converts "minutes after start of infusion" patterns to numeric hours.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "60 MIN AFTER START INF" → 1
#'
#' Converts minutes to hours by dividing by 60.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with MIN AFTER START patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_min_after_start <- function(xxtpt, result, na_idx) {
  min_after_start_pattern <- regex(
    "^(?<value>\\d+(?:\\.\\d+)?)\\s*m(?:in|inute)?s?\\s+after\\s+start\\s+inf$",
    ignore_case = TRUE,
    comments = TRUE
  )
  min_after_start_matches <- str_match(xxtpt, min_after_start_pattern)
  min_after_start_idx <- !is.na(min_after_start_matches[, 1]) &
    is.na(result) & !na_idx
  if (any(min_after_start_idx)) {
    result[min_after_start_idx] <- as.numeric(
      min_after_start_matches[min_after_start_idx, "value"]
    ) / 60
  }
  result
}

#' Convert MIN PRE EOI Patterns
#'
#' Converts "minutes pre end of infusion" patterns to negative numeric hours
#' (time before end of infusion).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "10MIN PRE EOI" → -0.1667
#'
#' Returns negative values to indicate time before end of infusion.
#' Converts minutes to hours by dividing by 60.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with MIN PRE EOI patterns converted to negative hours
#'
#' @keywords internal
#' @noRd
convert_min_pre_eoi <- function(xxtpt, result, na_idx) {
  min_pre_eoi_pattern <- regex(
    "^(?<value>\\d+(?:\\.\\d+)?)\\s*m(?:in|inute)?s?\\s+pre\\s+eoi$",
    ignore_case = TRUE,
    comments = TRUE
  )
  min_pre_eoi_matches <- str_match(xxtpt, min_pre_eoi_pattern)
  min_pre_eoi_idx <- !is.na(min_pre_eoi_matches[, 1]) & is.na(result) & !na_idx
  if (any(min_pre_eoi_idx)) {
    result[min_pre_eoi_idx] <- -as.numeric(
      min_pre_eoi_matches[min_pre_eoi_idx, "value"]
    ) / 60
  }
  result
}

#' Convert Simple Time Unit Patterns
#'
#' Converts simple hours-only and minutes-only patterns to numeric hours.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes and converts:
#' * Hours only: "1H", "2 hours", "4hr Post-dose" → 1, 2, 4
#' * Minutes only: "30M", "45 min", "30 MIN POST" → 0.5, 0.75, 0.5
#'
#' Processes hours first, then minutes. This function is called last to avoid
#' matching patterns that should be handled by more specific functions.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with simple unit patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_simple_units <- function(xxtpt, result, na_idx) {
  # Hours only
  hours_pattern <- regex(
    paste0(
      "^(?<hours>\\d+(?:\\.\\d+)?)\\s*h(?:r|our)?s?",
      "(?:\\s+post(?:\\s*-?\\s*dose)?)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hours_matches <- str_match(xxtpt, hours_pattern)
  hours_idx <- !is.na(hours_matches[, 1]) & is.na(result) & !na_idx
  if (any(hours_idx)) {
    result[hours_idx] <- as.numeric(hours_matches[hours_idx, "hours"])
  }

  # Minutes only
  minutes_pattern <- regex(
    paste0(
      "^(?<minutes>\\d+(?:\\.\\d+)?)\\s*m(?:in|inute)?s?",
      "(?:\\s+post(?:\\s*-?\\s*dose)?)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  minutes_matches <- str_match(xxtpt, minutes_pattern)
  minutes_idx <- !is.na(minutes_matches[, 1]) & is.na(result) & !na_idx
  if (any(minutes_idx)) {
    result[minutes_idx] <- as.numeric(
      minutes_matches[minutes_idx, "minutes"]
    ) / 60
  }

  result
}
