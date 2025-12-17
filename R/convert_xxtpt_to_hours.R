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
#' **Patterns Returning NA with Warning:**
#' * Time ranges with direction: `"0-4H PRIOR START OF INFUSION"`
#'
#' **Time-based Conversions:**
#' * **Days**: `"Day 1"` → 24, `"30 DAYS AFTER LAST"` → 720
#' * **Hours + Minutes**: `"1H30M"` → 1.5
#' * **Time Ranges**: `"0-6h"` → 6 (end value)
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
#' * Time range patterns with direction (with warning)
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
#' # Time ranges and combined units
#' convert_xxtpt_to_hours(c("1H30M", "0-6h Post-dose", "30 DAYS AFTER LAST"))
#'
#' # Custom infusion duration (2 hours)
#' # Note: "1 HOUR POST EOI" = 2 (infusion) + 1 (post) = 3 hours
#' convert_xxtpt_to_hours(
#'   c("EOI", "End of Infusion", "1 HOUR POST EOI"),
#'   infusion_duration = 2
#' )
#'
#' # After end patterns equal infusion_duration
#' convert_xxtpt_to_hours(c("AFTER END OF INFUSION", "AFTER END OF TREATMENT"))
#'
#' # Range patterns with direction return NA with warning
#' convert_xxtpt_to_hours(c("0-4H PRIOR START OF INFUSION"))
convert_xxtpt_to_hours <- function(xxtpt, infusion_duration = 1) {
  assert_character_vector(xxtpt)
  assert_numeric_vector(infusion_duration, length = 1)

  # Additional validation for positive value
  if (infusion_duration <= 0) {
    cli_abort(
      "{.arg infusion_duration} must be positive, but is {infusion_duration}."
    )
  }

  # Handle empty input
  if (length(xxtpt) == 0) {
    return(numeric(0))
  }

  # Trim whitespace from all inputs
  xxtpt <- trimws(xxtpt)

  # Initialize result vector with NA
  result <- rep(NA_real_, length(xxtpt))

  # Handle NA inputs (after trimming)
  na_idx <- is.na(xxtpt)

  # Track warnings for ambiguous/unsupported patterns
  range_values <- character(0)

  # Check special cases first (exact matches, case-insensitive) ----

  # Screening
  screening_pattern <- regex("^screening$", ignore_case = TRUE)
  screening_idx <- str_detect(xxtpt, screening_pattern) & !na_idx
  result[screening_idx] <- -1

  # Pre-dose, Predose, Pre-infusion, Pre-inf, Infusion, 0H, 0 H -> 0
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

  # Morning, Evening -> NA (already NA, no action needed)

  # Check for range patterns with direction (return NA with warning) ----
  # Only patterns like "0-4H PRIOR START OF INFUSION" should warn
  # Simple ranges like "0-6h Post-dose" are valid and handled below
  range_check_pattern <- regex(
    paste0(
      "^\\d+\\s*-\\s*\\d+\\s*h(?:r|our)?s?\\s+",
      "(?:prior|post)\\s+(?:start|end)"
    ),
    ignore_case = TRUE
  )
  range_check_idx <- str_detect(xxtpt, range_check_pattern) &
    is.na(result) & !na_idx
  if (any(range_check_idx)) {
    range_values <- c(range_values, unique(xxtpt[range_check_idx]))
  }

  # Check days (convert to hours by multiplying by 24) ----
  # Matches: "2D", "2 days", "Day 2", "2 day", "1.5 days", "Day 1"
  # Also: "30 DAYS AFTER LAST"
  days_pattern <- regex(
    paste0(
      # optional "day " prefix, then number
      "^(?:day\\s+)?(?<days>\\d+(?:\\.\\d+)?)\\s*",
      # OPTIONAL days suffix
      "(?:d|day|days)?",
      # optional "after last" or post-dose suffix
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

  # Check hours+minutes combinations ----
  # Matches: "1H30M", "1 hour 30 min", "1h 30m", "2HR15MIN"
  hm_pattern <- regex(
    paste0(
      # hours number
      "^(?<hours>\\d+(?:\\.\\d+)?)\\s*",
      # hours
      "h(?:r|our)?s?",
      # optional space
      "\\s*",
      # minutes number
      "(?<minutes>\\d+(?:\\.\\d+)?)\\s*",
      # minutes
      "m(?:in|inute)?s?",
      # optional post-dose suffix
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

  # Check time ranges (e.g., "0-6h Post-dose" -> 6, "0.5 - 6.5h" -> 6.5) ----
  # Process before simple hours to avoid conflicts
  range_pattern <- regex(
    paste0(
      # range with end value captured, spaces allowed
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      # hours
      "h(?:r|our)?s?",
      # optional post-dose suffix
      "(?:\\s+post(?:\\s*-?\\s*dose)?)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  range_matches <- str_match(xxtpt, range_pattern)
  range_idx <- !is.na(range_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_idx)) {
    result[range_idx] <- as.numeric(range_matches[range_idx, "end"])
  }

  # Check "X MIN/HOUR PREDOSE" (negative time) ----
  predose_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "predose$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  predose_matches <- str_match(xxtpt, predose_pattern)
  predose_idx <- !is.na(predose_matches[, 1]) & is.na(result) & !na_idx
  if (any(predose_idx)) {
    time_value <- as.numeric(predose_matches[predose_idx, "value"])
    unit <- tolower(predose_matches[predose_idx, "unit"])
    # Check if unit starts with 'm' (minutes)
    is_minutes <- substr(unit, 1, 1) == "m"
    result[predose_idx] <- if_else(
      is_minutes,
      -time_value / 60,
      -time_value
    )
  }

  # Check "X HOUR/MIN POST EOI/EOT" patterns ----
  post_eoi_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "(?<unit>m(?:in|inute)?|h(?:r|our)?)s?\\s+",
      "post\\s+(?:eoi|eot)$"
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
    # Add infusion_duration since this is time AFTER end of infusion
    result[post_eoi_idx] <- if_else(
      is_minutes,
      infusion_duration + time_value / 60,
      infusion_duration + time_value
    )
  }

  # Check "XMIN/XHOUR AFTER END OF INFUSION/TREATMENT" (compact format) ----
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
    # Add infusion_duration since this is time AFTER end of infusion
    result[after_end_idx] <- if_else(
      is_minutes,
      infusion_duration + time_value / 60,
      infusion_duration + time_value
    )
  }

  # Check "X HR POST INF/EOI/EOT" patterns ----
  hr_post_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "hr\\s+post\\s+(?:inf|eoi|eot)$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hr_post_matches <- str_match(xxtpt, hr_post_pattern)
  hr_post_idx <- !is.na(hr_post_matches[, 1]) & is.na(result) & !na_idx
  if (any(hr_post_idx)) {
    # Add infusion_duration since this is time AFTER end of infusion
    result[hr_post_idx] <- infusion_duration +
      as.numeric(hr_post_matches[hr_post_idx, "value"])
  }

  # Check "XH PRIOR/POST START OF INFUSION" (compact H) ----
  h_start_inf_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)h\\s+",
      "(?<direction>prior|post)\\s+",
      "start\\s+of\\s+infusion$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  h_start_inf_matches <- str_match(xxtpt, h_start_inf_pattern)
  h_start_inf_idx <- !is.na(h_start_inf_matches[, 1]) &
    is.na(result) & !na_idx
  if (any(h_start_inf_idx)) {
    time_value <- as.numeric(h_start_inf_matches[h_start_inf_idx, "value"])
    direction <- tolower(h_start_inf_matches[h_start_inf_idx, "direction"])
    result[h_start_inf_idx] <- if_else(
      direction == "prior",
      -time_value,
      time_value
    )
  }

  # Check "X MIN AFTER START INF" ----
  min_after_start_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "m(?:in|inute)?s?\\s+",
      "after\\s+start\\s+inf$"
    ),
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

  # Check "XMIN PRE EOI" (compact, pre end of infusion) ----
  min_pre_eoi_pattern <- regex(
    paste0(
      "^(?<value>\\d+(?:\\.\\d+)?)\\s*",
      "m(?:in|inute)?s?\\s+",
      "pre\\s+eoi$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  min_pre_eoi_matches <- str_match(xxtpt, min_pre_eoi_pattern)
  min_pre_eoi_idx <- !is.na(min_pre_eoi_matches[, 1]) &
    is.na(result) & !na_idx
  if (any(min_pre_eoi_idx)) {
    result[min_pre_eoi_idx] <- -as.numeric(
      min_pre_eoi_matches[min_pre_eoi_idx, "value"]
    ) / 60
  }

  # Check hours only ----
  # Matches: "1H", "2 hours", "0.5HR", "1h Post-dose", "1H Post dose",
  # "8 hour", "12 HOURS". Also: "1 HOUR POST"
  hours_pattern <- regex(
    paste0(
      # number
      "^(?<hours>\\d+(?:\\.\\d+)?)\\s*",
      # hours
      "h(?:r|our)?s?",
      # optional post/post-dose suffix (flexible spacing)
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

  # Check minutes only ----
  # Matches: "30M", "45 min", "30 Min Post-dose", "30 MIN Post dose",
  # "2.5 Min". Also: "X MIN POST"
  minutes_pattern <- regex(
    paste0(
      # number
      "^(?<minutes>\\d+(?:\\.\\d+)?)\\s*",
      # minutes
      "m(?:in|inute)?s?",
      # optional post/post-dose suffix (flexible spacing)
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

  # Issue warnings for ambiguous patterns ----
  if (length(range_values) > 0) {
    cli_warn(
      c(
        paste(
          "Time range patterns with direction cannot be converted",
          "to single numeric hours."
        ),
        "i" = "Returning NA for: {.val {range_values}}"
      )
    )
  }

  result
}
