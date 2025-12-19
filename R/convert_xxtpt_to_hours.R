#' Convert `XXTPT` Strings to Hours
#'
#' Converts CDISC timepoint strings (e.g., `PCTPT`, `VSTPT`, `EGTPT`,
#' `ISTPT`, `LBTPT`) into numeric hours for analysis. The function handles
#' common dose-centric formats including pre-dose, post-dose (hours/minutes),
#' days, time ranges, and treatment-related time markers.
#'
#' @param xxtpt A character vector of timepoint descriptions from SDTM
#'   `--TPT` variables (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`).
#'   Can contain `NA` values.
#'
#' @permitted character vector
#'
#' @param treatment_duration Numeric value(s) specifying the duration of treatment
#'   in hours. Used to convert "EOI/EOT" (End of Infusion/Treatment) patterns
#'   and patterns describing time after end of treatment. Must be non-negative.
#'   Can be either:
#'   * A single value (used for all timepoints), or
#'   * A vector of the same length as `xxtpt` (one value per timepoint)
#'
#'   Default is 0 hours (for instantaneous treatments like oral medications).
#'
#' @permitted numeric scalar or vector (non-negative)
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
#' * `"Pre-dose"`, `"Predose"`, `"Pre-treatment"`, `"Pre-infusion"`,
#'   `"Pre-inf"`, `"Before"`, `"Infusion"`, `"0H"` → 0
#' * `"EOI"`, `"EOT"`, `"End of Infusion"`, `"End of Treatment"`,
#'   `"After End of Infusion"`, `"After End of Treatment"` → `treatment_duration`
#'   (default: 0)
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
#' * `"0-4H AFTER EOI"` with midpoint and treatment_duration=1 → 3 (1 + 2)
#' * `"0-4H EOT"` with midpoint and treatment_duration=0 → 2
#' * `"4-8H AFTER END OF INFUSION"` with midpoint and treatment_duration=1 → 7 (1 + 6)
#'
#' **Time-based Conversions:**
#' * **Days**: `"Day 1"` → 24, `"30 DAYS AFTER LAST"` → 720
#' * **Hours + Minutes**: `"1H30M"` → 1.5
#' * **Hours**: `"2 hours"` → 2, `"1 HOUR POST"` → 1
#' * **Minutes**: `"30M"` → 0.5, `"30 MIN POST"` → 0.5
#' * **Predose**: `"5 MIN PREDOSE"` → -0.0833
#' * **Before treatment**: `"5 MIN BEFORE"` → -0.0833
#' * **Post EOI/EOT**: `"1 HOUR POST EOI"` → treatment_duration + 1,
#'   `"24 HR POST INF"` → treatment_duration + 24,
#'   `"1 HOUR AFTER EOT"` → treatment_duration + 1
#' * **After end**: `"30MIN AFTER END OF INFUSION"` → treatment_duration + 0.5
#' * **Start of infusion/treatment**: `"8H PRIOR START OF INFUSION"` → -8,
#'   `"8H BEFORE START OF TREATMENT"` → -8
#' * **Pre EOI/EOT**: `"10MIN PRE EOI"` → treatment_duration - 1/6,
#'   `"10MIN BEFORE EOT"` → treatment_duration - 1/6
#'
#' **Supported Unit Formats:**
#' * Hours: H, h, HR, hr, HOUR, hour (with optional plurals)
#' * Minutes: M, m, MIN, min, MINUTE, minute (with optional plurals)
#' * Days: D, d, DAY, day (with optional plurals)
#' * Flexible whitespace and optional "Post-dose", "POST", "After last"
#'   suffixes
#'
#' **Understanding POST/AFTER Patterns:**
#'
#' It's important to distinguish between patterns relative to treatment **start**
#' versus treatment **end**:
#'
#' * **Relative to START** (treatment_duration NOT added):
#'   - `"1H POST"`, `"1H AFTER"`, `"30M POST"` → Time from dose/treatment start
#'   - These patterns assume treatment starts at time 0
#'   - Example: `"1H POST"` → 1 hour (regardless of treatment_duration)
#'
#' * **Relative to END** (treatment_duration IS added):
#'   - `"1H POST EOI"`, `"1H AFTER EOT"`, `"1H POST INFUSION"` → Time from
#'     treatment end
#'   - These patterns account for when treatment ends (start + duration)
#'   - Example: `"1H POST EOI"` with treatment_duration=2 → 3 hours (2 + 1)
#'
#' This distinction follows standard pharmacokinetic conventions where "post-dose"
#' refers to time from treatment initiation, while "post end of infusion" refers
#' to time from treatment completion.
#'
#' **Vectorized Treatment Duration:**
#'
#' When `treatment_duration` is a vector, each timepoint uses its corresponding
#' treatment duration value. This is useful when different records have different
#' treatment durations (e.g., different infusion lengths).
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
#'   "Pre-treatment",
#'   "Before",
#'   "30M",
#'   "1H",
#'   "2H POSTDOSE",
#'   "Day 1"
#' ))
#'
#' # Predose/before patterns (negative times)
#' convert_xxtpt_to_hours(c("5 MIN PREDOSE", "1 HOUR BEFORE"))
#'
#' # Treatment-related patterns with default treatment_duration = 0 (oral)
#' convert_xxtpt_to_hours(c(
#'   "EOT",
#'   "1 HOUR POST EOT",
#'   "1 HOUR AFTER EOT",
#'   "After End of Treatment"
#' ))
#'
#' # Infusion-related patterns with treatment_duration = 1 hour (scalar)
#' convert_xxtpt_to_hours(
#'   c(
#'     "EOI",
#'     "1 HOUR POST EOI",
#'     "24 HR POST INF",
#'     "30MIN AFTER END OF INFUSION",
#'     "8H PRIOR START OF INFUSION",
#'     "10MIN PRE EOI"
#'   ),
#'   treatment_duration = 1
#' )
#'
#' # Vectorized treatment_duration - different durations per timepoint
#' convert_xxtpt_to_hours(
#'   c("EOI", "1 HOUR POST EOI", "EOI", "1 HOUR POST EOI"),
#'   treatment_duration = c(1, 1, 2, 2)
#' )
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
#' # Demonstrating POST vs POST EOI distinction
#' # POST alone = relative to treatment START (no duration added)
#' convert_xxtpt_to_hours(
#'   c("Pre-dose", "1H POST", "2H POST", "4H POST"),
#'   treatment_duration = 2 # 2-hour infusion
#' )
#' # Returns: 0, 1, 2, 4 (treatment_duration NOT added)
#'
#' # POST EOI = relative to treatment END (duration added)
#' convert_xxtpt_to_hours(
#'   c("Pre-dose", "EOI", "1H POST EOI", "2H POST EOI"),
#'   treatment_duration = 2 # 2-hour infusion
#' )
#' # Returns: 0, 2, 3, 4 (treatment_duration IS added)
#'
#' # Comparing the two patterns side by side
#' convert_xxtpt_to_hours(
#'   c("1H POST", "1H POST EOI", "1H POST INFUSION"),
#'   treatment_duration = 2
#' )
#' # Returns: 1, 3, 3
#' # "1H POST" = 1 hour from start
#' # "1H POST EOI" = 1 hour after 2-hour infusion ends = 3 hours from start
#' # "1H POST INFUSION" = same as POST EOI
convert_xxtpt_to_hours <- function(xxtpt,
                                   treatment_duration = 0,
                                   range_method = "midpoint") {
  # Validate inputs
  assert_character_vector(xxtpt)
  assert_numeric_vector(treatment_duration)
  assert_character_vector(range_method, values = c("start", "end", "midpoint"))

  # Validate treatment_duration length
  if (length(treatment_duration) != 1 && length(treatment_duration) != length(xxtpt)) {
    cli_abort(
      c(
        "{.arg treatment_duration} must be either:",
        "i" = "A single value (used for all timepoints), or",
        "i" = "A vector of length {length(xxtpt)} (one value per timepoint)",
        "x" = "You've supplied a vector of length {length(treatment_duration)}
        for {length(xxtpt)} timepoint{?s}."
      )
    )
  }

  # Recycle if scalar
  if (length(treatment_duration) == 1 && length(xxtpt) > 1) {
    treatment_duration <- rep(treatment_duration, length(xxtpt))
  }

  # Validate all values are non-negative
  if (any(treatment_duration < 0, na.rm = TRUE)) {
    cli_abort(
      "{.arg treatment_duration} must be non-negative for all values."
    )
  }

  if (length(xxtpt) == 0) {
    return(numeric(0))
  }

  xxtpt <- trimws(xxtpt)
  result <- rep(NA_real_, length(xxtpt))
  na_idx <- is.na(xxtpt)

  # Process patterns in order of specificity
  result <- convert_special_cases(xxtpt, result, na_idx, treatment_duration)
  result <- convert_time_units(xxtpt, result, na_idx)
  result <- convert_ranges_eot(xxtpt, result, na_idx, treatment_duration, range_method)
  result <- convert_ranges_after_end(xxtpt, result, na_idx, treatment_duration, range_method)
  result <- convert_ranges(xxtpt, result, na_idx, range_method)
  result <- convert_treatment_patterns(xxtpt, result, na_idx, treatment_duration)
  result <- convert_simple_units(xxtpt, result, na_idx)

  result
}

# Helper Functions ----

#' Convert Special Case Patterns
#'
#' Converts special case timepoint patterns to numeric hours including screening,
#' pre-dose, pre-treatment, before, and end of treatment markers.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes and converts the following patterns:
#' * "Screening" → -1
#' * "Pre-dose", "Predose", "Pre-treatment", "Pre-infusion", "Pre-inf",
#'   "Before", "Infusion", "0H" → 0
#' * "EOI", "EOT", "End of Infusion", "End of Treatment",
#'   "After End of Infusion", "After End of Treatment" → treatment_duration
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with special case patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_special_cases <- function(xxtpt, result, na_idx, treatment_duration) {
  # Screening
  screening_pattern <- regex("^screening$", ignore_case = TRUE)
  screening_idx <- str_detect(xxtpt, screening_pattern) & !na_idx
  result[screening_idx] <- -1

  # Pre-dose, Predose, Pre-treatment, Pre-infusion, Pre-inf, Before, Infusion, 0H -> 0
  zero_pattern <- regex(
    paste0(
      "^(pre-?(?:dose|treatment|inf(?:usion)?)|",
      "before|infusion|0\\s*h(?:r|our)?s?)$"
    ),
    ignore_case = TRUE
  )
  zero_idx <- str_detect(xxtpt, zero_pattern) & is.na(result) & !na_idx
  result[zero_idx] <- 0

  # EOI, EOT, End of Infusion/Treatment, After End of Infusion/Treatment -> treatment_duration
  eot_pattern <- regex(
    paste0(
      "^(eo[it]|end\\s+of\\s+(?:infusion|treatment)|",
      "after\\s+end\\s+of\\s+(?:infusion|treatment))$"
    ),
    ignore_case = TRUE
  )
  eot_idx <- str_detect(xxtpt, eot_pattern) & is.na(result) & !na_idx
  result[eot_idx] <- treatment_duration[eot_idx]

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
      "(?:\\s+(?:post|after)(?:\\s*-?\\s*dose)?)?$"
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
#' * Ranges with direction: "0-4H PRIOR START OF INFUSION" → -2 (with midpoint),
#'   "0-4H BEFORE START OF TREATMENT" → -2 (with midpoint)
#' * Simple ranges: "0-6h Post-dose" → 3 (with midpoint)
#'
#' Direction affects sign: "PRIOR/BEFORE" results in negative values, "POST/AFTER"
#' in positive. Accepts both "INFUSION" and "TREATMENT" terminology.
#'
#' Processes ranges with direction before simple ranges to ensure correct matching.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with range patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_ranges <- function(xxtpt, result, na_idx, range_method) {
  # Ranges with direction (PRIOR/BEFORE/POST/AFTER START/END)
  range_dir_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?\\s+",
      "(?<direction>prior|before|post|after)\\s+",
      "(?:start|end)(?:\\s+of\\s+(?:infusion|treatment))?"
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
      direction %in% c("prior", "before"),
      -range_val,
      range_val
    )
  }

  # Simple time ranges
  range_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?(?:\\s+(?:post|after)(?:\\s*-?\\s*dose)?)?$"
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

#' Convert Range Patterns Relative to EOI/EOT
#'
#' Converts time range patterns relative to end of infusion/treatment to numeric
#' hours using specified range method, adding the treatment duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#' @param range_method Method for selecting value from range: "start", "midpoint",
#'   or "end"
#'
#' @details
#' Recognizes and converts the following range patterns:
#' * "0-4H AFTER EOI" → treatment_duration + range_value
#' * "0-4H EOI" → treatment_duration + range_value
#' * "0-4H AFTER EOT" → treatment_duration + range_value
#' * "0-4H EOT" → treatment_duration + range_value
#'
#' With midpoint method, "0-4H AFTER EOI" with treatment_duration=1 → 1 + 2 = 3
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with EOI/EOT range patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_ranges_eot <- function(xxtpt, result, na_idx, treatment_duration, range_method) {
  # Ranges relative to EOI/EOT (with optional AFTER)
  range_eot_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?\\s+",
      "(?:after\\s+)?eo[it]$"
    ),
    ignore_case = TRUE
  )
  range_eot_matches <- str_match(xxtpt, range_eot_pattern)
  range_eot_idx <- !is.na(range_eot_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_eot_idx)) {
    start_val <- as.numeric(range_eot_matches[range_eot_idx, "start"])
    end_val <- as.numeric(range_eot_matches[range_eot_idx, "end"])
    range_val <- calculate_range_value(start_val, end_val, range_method)
    result[range_eot_idx] <- treatment_duration[range_eot_idx] + range_val
  }

  result
}

#' Convert Range Patterns After End of Infusion/Treatment
#'
#' Converts time range patterns after end of infusion/treatment (long form) to
#' numeric hours using specified range method, adding the treatment duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#' @param range_method Method for selecting value from range: "start", "midpoint",
#'   or "end"
#'
#' @details
#' Recognizes and converts the following range patterns:
#' * "4-8H AFTER END OF INFUSION" → treatment_duration + range_value
#' * "4-8H AFTER END OF TREATMENT" → treatment_duration + range_value
#'
#' With midpoint method, "4-8H AFTER END OF INFUSION" with treatment_duration=1
#' → 1 + 6 = 7
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with after end range patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_ranges_after_end <- function(xxtpt, result, na_idx, treatment_duration, range_method) {
  # Ranges after end of infusion/treatment (long form)
  range_after_end_pattern <- regex(
    paste0(
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      "h(?:r|our)?s?\\s+",
      "after\\s+end\\s+of\\s+(?:infusion|treatment)$"
    ),
    ignore_case = TRUE
  )
  range_after_end_matches <- str_match(xxtpt, range_after_end_pattern)
  range_after_end_idx <- !is.na(range_after_end_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_after_end_idx)) {
    start_val <- as.numeric(range_after_end_matches[range_after_end_idx, "start"])
    end_val <- as.numeric(range_after_end_matches[range_after_end_idx, "end"])
    range_val <- calculate_range_value(start_val, end_val, range_method)
    result[range_after_end_idx] <- treatment_duration[range_after_end_idx] + range_val
  }

  result
}

#' Convert Treatment-Related Patterns
#'
#' Orchestrates conversion of all treatment-related timepoint patterns by calling
#' specialized helper functions for each pattern type.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Calls helper functions in sequence to convert:
#' * Pre-treatment/predose patterns (negative times)
#' * Post EOI/EOT patterns
#' * After end of infusion/treatment patterns
#' * HR POST INF patterns
#' * POST INFUSION/INF patterns
#' * Start of infusion/treatment patterns
#' * MIN AFTER START patterns
#' * MIN PRE/BEFORE EOI/EOT patterns
#'
#' @return Updated numeric vector with all treatment-related patterns converted
#'
#' @keywords internal
#' @noRd
convert_treatment_patterns <- function(xxtpt, result, na_idx, treatment_duration) {
  result <- convert_predose_patterns(xxtpt, result, na_idx)
  result <- convert_post_eot_patterns(xxtpt, result, na_idx, treatment_duration)
  result <- convert_after_end_patterns(xxtpt, result, na_idx, treatment_duration)
  result <- convert_hr_post_patterns(xxtpt, result, na_idx, treatment_duration)
  result <- convert_post_inf_patterns(xxtpt, result, na_idx, treatment_duration)
  result <- convert_start_patterns(xxtpt, result, na_idx)
  result <- convert_min_after_start(xxtpt, result, na_idx)
  result <- convert_min_pre_eot(xxtpt, result, na_idx, treatment_duration)
  result
}

#' Convert Predose/Before Patterns
#'
#' Converts predose and before timepoint patterns to negative numeric hours
#' (time before dose/treatment).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "5 MIN PREDOSE" → -0.0833 (negative 5 minutes)
#' * "1 HOUR BEFORE" → -1 (negative 1 hour)
#'
#' Returns negative values to indicate time before dose/treatment.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with predose/before patterns converted to negative hours
#'
#' @keywords internal
#' @noRd
convert_predose_patterns <- function(xxtpt, result, na_idx) {
  predose_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "(?:predose|before)$"
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
#' Converts "post/after end of infusion/treatment" patterns to numeric hours,
#' adding the treatment duration to account for time after treatment completion.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes patterns like:
#' * "1 HOUR POST EOI" → treatment_duration + 1
#' * "30 MIN AFTER EOT" → treatment_duration + 0.5
#' * "1 HOUR POST EOT" → treatment_duration + 1
#'
#' Adds treatment_duration because these patterns describe time AFTER the end of
#' treatment, so total time from treatment start includes the treatment duration
#' plus the additional time.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with post EOI/EOT patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_post_eot_patterns <- function(xxtpt, result, na_idx, treatment_duration) {
  post_eot_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "(?:post|after)\\s+eo[it]$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  post_eot_matches <- str_match(xxtpt, post_eot_pattern)
  post_eot_idx <- !is.na(post_eot_matches[, 1]) & is.na(result) & !na_idx
  if (any(post_eot_idx)) {
    time_value <- as.numeric(post_eot_matches[post_eot_idx, "value"])
    unit <- tolower(post_eot_matches[post_eot_idx, "unit"])
    is_minutes <- substr(unit, 1, 1) == "m"
    result[post_eot_idx] <- if_else(
      is_minutes,
      treatment_duration[post_eot_idx] + time_value / 60,
      treatment_duration[post_eot_idx] + time_value
    )
  }
  result
}

#' Convert After End of Infusion/Treatment Patterns
#'
#' Converts "after end of infusion/treatment" patterns (compact format) to numeric
#' hours, adding the treatment duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes patterns like:
#' * "30MIN AFTER END OF INFUSION" → treatment_duration + 0.5
#' * "1 HOUR AFTER END OF TREATMENT" → treatment_duration + 1
#'
#' Adds treatment_duration because these patterns describe time AFTER the end of
#' treatment.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with after end patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_after_end_patterns <- function(xxtpt, result, na_idx, treatment_duration) {
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
      treatment_duration[after_end_idx] + time_value / 60,
      treatment_duration[after_end_idx] + time_value
    )
  }
  result
}

#' Convert HR POST INF/EOI/EOT Patterns
#'
#' Converts "HR POST INF/EOI/EOT" patterns to numeric hours, adding the treatment
#' duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes patterns like:
#' * "24 HR POST INF" → treatment_duration + 24
#' * "24 HR POST EOI" → treatment_duration + 24
#' * "24 HR POST EOT" → treatment_duration + 24
#'
#' Adds treatment_duration because these patterns describe time AFTER the end of
#' treatment.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with HR POST patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_hr_post_patterns <- function(xxtpt, result, na_idx, treatment_duration) {
  hr_post_pattern <- regex(
    "^(?<value>\\d+(?:\\.\\d+)?)\\s*hr\\s+post\\s+(?:inf|eo[it])$",
    ignore_case = TRUE,
    comments = TRUE
  )
  hr_post_matches <- str_match(xxtpt, hr_post_pattern)
  hr_post_idx <- !is.na(hr_post_matches[, 1]) & is.na(result) & !na_idx
  if (any(hr_post_idx)) {
    result[hr_post_idx] <- treatment_duration[hr_post_idx] +
      as.numeric(hr_post_matches[hr_post_idx, "value"])
  }
  result
}

#' Convert POST INFUSION/INF Patterns
#'
#' Converts "post infusion/inf" patterns (without EOI/EOT) to numeric hours,
#' adding the treatment duration.
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes patterns like:
#' * "30 MIN POST INFUSION" → treatment_duration + 0.5
#' * "1 HOUR POST INF" → treatment_duration + 1
#'
#' Adds treatment_duration because these patterns describe time AFTER the end of
#' infusion/treatment.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with post infusion patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_post_inf_patterns <- function(xxtpt, result, na_idx, treatment_duration) {
  post_inf_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "post\\s+inf(?:usion)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  post_inf_matches <- str_match(xxtpt, post_inf_pattern)
  post_inf_idx <- !is.na(post_inf_matches[, 1]) & is.na(result) & !na_idx
  if (any(post_inf_idx)) {
    time_value <- as.numeric(post_inf_matches[post_inf_idx, "value"])
    unit <- tolower(post_inf_matches[post_inf_idx, "unit"])
    is_minutes <- substr(unit, 1, 1) == "m"
    result[post_inf_idx] <- if_else(
      is_minutes,
      treatment_duration[post_inf_idx] + time_value / 60,
      treatment_duration[post_inf_idx] + time_value
    )
  }
  result
}

#' Convert Start of Infusion/Treatment Patterns
#'
#' Converts patterns relative to the start of infusion/treatment to numeric hours,
#' with sign based on direction (prior/before = negative, post/after = positive).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#'
#' @details
#' Recognizes patterns like:
#' * "8H PRIOR START OF INFUSION" → -8
#' * "8H BEFORE START OF TREATMENT" → -8
#' * "8H POST START OF INFUSION" → 8
#' * "8H AFTER START OF TREATMENT" → 8
#'
#' Direction affects sign: "PRIOR/BEFORE" results in negative values,
#' "POST/AFTER" in positive. Accepts both "INFUSION" and "TREATMENT" terminology.
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with start of treatment patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_start_patterns <- function(xxtpt, result, na_idx) {
  h_start_treatment_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)h\\s+",
      "(?<direction>prior|before|post|after)\\s+",
      "start\\s+of\\s+(?:infusion|treatment)$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  h_start_treatment_matches <- str_match(xxtpt, h_start_treatment_pattern)
  h_start_treatment_idx <- !is.na(h_start_treatment_matches[, 1]) &
    is.na(result) & !na_idx
  if (any(h_start_treatment_idx)) {
    time_value <- as.numeric(h_start_treatment_matches[h_start_treatment_idx, "value"])
    direction <- tolower(h_start_treatment_matches[h_start_treatment_idx, "direction"])
    result[h_start_treatment_idx] <- if_else(
      direction %in% c("prior", "before"),
      -time_value,
      time_value
    )
  }
  result
}

#' Convert MIN AFTER START INF Patterns
#'
#' Converts "minutes after start of infusion/treatment" patterns to numeric hours.
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

#' Convert MIN PRE/BEFORE EOI/EOT Patterns
#'
#' Converts "minutes pre/before end of infusion/treatment" patterns to hours
#' relative to the end of treatment (treatment_duration minus the time).
#'
#' @param xxtpt Character vector of timepoint descriptions (trimmed, no leading/
#'   trailing whitespace)
#' @param result Numeric vector of results (partially filled, may contain NA)
#' @param na_idx Logical vector indicating which positions in xxtpt are NA
#' @param treatment_duration Duration of treatment in hours (non-negative numeric
#'   vector, same length as xxtpt)
#'
#' @details
#' Recognizes patterns like:
#' * "10MIN PRE EOI" → treatment_duration - 1/6
#' * "10MIN BEFORE EOT" → treatment_duration - 1/6
#'
#' Returns treatment_duration minus the specified time to indicate time before
#' end of treatment. For example, with treatment_duration = 2 hours, "10MIN PRE EOI"
#' returns 2 - 10/60 = 1.833 hours (10 minutes before the 2-hour infusion ends).
#'
#' Only updates result for positions where result is currently NA and xxtpt is not NA.
#'
#' @return Updated numeric vector with MIN PRE/BEFORE EOI/EOT patterns converted to hours
#'
#' @keywords internal
#' @noRd
convert_min_pre_eot <- function(xxtpt, result, na_idx, treatment_duration) {
  min_pre_eot_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "m(?:in|inute)?s?\\s+",
      "(?:pre|before)\\s+eo[it]$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  min_pre_eot_matches <- str_match(xxtpt, min_pre_eot_pattern)
  min_pre_eot_idx <- !is.na(min_pre_eot_matches[, 1]) & is.na(result) & !na_idx
  if (any(min_pre_eot_idx)) {
    result[min_pre_eot_idx] <- treatment_duration[min_pre_eot_idx] - as.numeric(
      min_pre_eot_matches[min_pre_eot_idx, "value"]
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
#' * Hours only: "1H", "2 hours", "4hr Post-dose", "4hr After-dose" → 1, 2, 4, 4
#' * Minutes only: "30M", "45 min", "30 MIN POST", "30 MIN AFTER" → 0.5, 0.75, 0.5, 0.5
#'
#' Processes hours first, then minutes. This function is called last to avoid
#' matching patterns that should be handled by more specific functions.
#'
#' Accepts both "POST" and "AFTER" terminology.
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
      "(?:\\s+(?:post|after)(?:\\s*-?\\s*dose)?)?$"
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
      "(?:\\s+(?:post|after)(?:\\s*-?\\s*dose)?)?$"
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
