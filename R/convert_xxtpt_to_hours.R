#' Convert `XXTPT` Strings to Hours
#'
#' Converts CDISC timepoint strings (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`)
#' into numeric hours for analysis. The function handles common dose-centric formats
#' including pre-dose, post-dose (hours/minutes), days, time ranges, and special
#' time markers.
#'
#' @param xxtpt A character vector of timepoint descriptions from SDTM `--TPT` variables
#'   (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`). Can contain `NA` values.
#'
#' @details
#' The function recognizes the following patterns (all case-insensitive):
#'
#' **Special Cases:**
#' * `"Screening"` → -1
#' * `"Pre-dose"`, `"Predose"`, `"Pre-infusion"`, `"Infusion"`, `"0H"`, `"0 H"` → 0
#' * `"EOI"`, `"End of Infusion"` → 1
#' * `"Morning"`, `"Evening"` → `NA_real_`
#' * Unrecognized values → `NA_real_`
#'
#' **Time-based Conversions:**
#' * **Days**: `"Day 1"`, `"2D"`, `"2 days"`, `"1.5 days"` → multiply by 24 (e.g., 24, 48, 36)
#' * **Hours + Minutes**: `"1H30M"`, `"1 hour 30 min"`, `"2HR15MIN"` → hours +
#' minutes/60 (e.g., 1.5, 2.25)
#' * **Time Ranges**: `"0-6h"`, `"6-12h Post-dose"`, `"0.5 - 6.5h"` → end value (e.g., 6, 12, 6.5)
#' * **Hours only**: `"1H"`, `"2 hours"`, `"0.5HR"`, `"4hr Post-dose"` → as-is (e.g., 1, 2, 0.5, 4)
#' * **Minutes only**: `"30M"`, `"45 min"`, `"2.5 Min"` → divide by 60 (e.g., 0.5, 0.75, 0.042)
#'
#' **Supported Unit Formats:**
#' * Hours: H, h, HR, hr, HOUR, hour, HOURS, hours
#' * Minutes: M, m, MIN, min, MINUTE, minute, MINUTES, minutes
#' * Days: D, d, DAY, day, DAYS, days
#' * All with optional plurals and flexible whitespace
#' * Optional "Post-dose" or "Post dose" suffix supported for all patterns
#'
#' @return A numeric vector of timepoints in hours. Returns `NA_real_` for:
#' * Input `NA` values
#' * Unrecognized timepoint formats
#' * Non-time descriptors (e.g., "Morning", "Evening")
#'
#' Returns `numeric(0)` for empty input.
#'
#' @keywords utils_fmt
#' @family utils_fmt
#'
#' @export
#'
#' @examples
#' convert_xxtpt_to_hours(c(
#'   "Screening", "Pre-dose", "EOI", "5 Min Post-dose", "10 MIN POSTDOSE", "30M",
#'   "1H30M", "2H POSTDOSE", "0-6h Post-dose", "Day 1", "24h", "48 HR"
#' ))
convert_xxtpt_to_hours <- function(xxtpt) {
  assert_character_vector(xxtpt)

  # Handle empty input
  if (length(xxtpt) == 0) {
    return(numeric(0))
  }

  # Initialize result vector with NA
  result <- rep(NA_real_, length(xxtpt))

  # Handle NA inputs
  na_idx <- is.na(xxtpt)

  # 1. Check special cases first (exact matches, case-insensitive)

  # Screening
  screening_pattern <- regex("^screening$", ignore_case = TRUE)
  screening_idx <- str_detect(xxtpt, screening_pattern) & !na_idx
  result[screening_idx] <- -1

  # Pre-dose, Predose, Pre-infusion, Infusion, 0H, 0 H -> 0
  zero_pattern <- regex(
    "^(pre-?dose|pre-?infusion|infusion|0\\s*h(?:r|our)?s?)$",
    ignore_case = TRUE
  )
  zero_idx <- str_detect(xxtpt, zero_pattern) & is.na(result) & !na_idx
  result[zero_idx] <- 0

  # EOI, End of Infusion -> 1
  eoi_pattern <- regex("^(eoi|end\\s+of\\s+infusion)$", ignore_case = TRUE)
  eoi_idx <- str_detect(xxtpt, eoi_pattern) & is.na(result) & !na_idx
  result[eoi_idx] <- 1

  # Morning, Evening -> NA (already NA, no action needed)

  # 2. Check days (convert to hours by multiplying by 24)
  # Matches: "2D", "2 days", "Day 2", "2 day", "1.5 days", "Day 1"
  days_pattern <- regex(
    paste0(
      # optional "day " prefix, then number
      "^(?:day\\s+)?(?<days>\\d+(?:\\.\\d+)?)\\s*",
      # OPTIONAL days suffix
      "(?:d|day|days)?",
      # optional post-dose suffix
      "(?:\\s+post-?dose)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  days_matches <- str_match(xxtpt, days_pattern)
  days_idx <- !is.na(days_matches[, 1]) & is.na(result) & !na_idx
  if (any(days_idx)) {
    result[days_idx] <- as.numeric(days_matches[days_idx, "days"]) * 24
  }

  # 3. Check hours+minutes combinations
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
      "(?:\\s+post-?dose)?$"
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

  # 4. Check time ranges (e.g., "0-6h Post-dose" -> 6, "0.5 - 6.5h" -> 6.5)
  # Process before simple hours to avoid conflicts
  range_pattern <- regex(
    paste0(
      # range with end value captured, spaces allowed
      "^(?<start>\\d+(?:\\.\\d+)?)\\s*-\\s*(?<end>\\d+(?:\\.\\d+)?)\\s*",
      # hours
      "h(?:r|our)?s?",
      # optional post-dose suffix
      "(?:\\s+post-?dose)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  range_matches <- str_match(xxtpt, range_pattern)
  range_idx <- !is.na(range_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_idx)) {
    result[range_idx] <- as.numeric(range_matches[range_idx, "end"])
  }

  # 5. Check hours only
  # Matches: "1H", "2 hours", "0.5HR", "1h Post-dose", "8 hour", "12 HOURS"
  hours_pattern <- regex(
    paste0(
      # number
      "^(?<hours>\\d+(?:\\.\\d+)?)\\s*",
      # hours
      "h(?:r|our)?s?",
      # optional post-dose suffix
      "(?:\\s+post-?dose)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hours_matches <- str_match(xxtpt, hours_pattern)
  hours_idx <- !is.na(hours_matches[, 1]) & is.na(result) & !na_idx
  if (any(hours_idx)) {
    result[hours_idx] <- as.numeric(hours_matches[hours_idx, "hours"])
  }

  # 6. Check minutes only
  # Matches: "30M", "45 min", "30 Min Post-dose", "2.5 Min"
  minutes_pattern <- regex(
    paste0(
      # number
      "^(?<minutes>\\d+(?:\\.\\d+)?)\\s*",
      # minutes
      "m(?:in|inute)?s?",
      # optional post-dose suffix
      "(?:\\s+post-?dose)?$"
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  minutes_matches <- str_match(xxtpt, minutes_pattern)
  minutes_idx <- !is.na(minutes_matches[, 1]) & is.na(result) & !na_idx
  if (any(minutes_idx)) {
    result[minutes_idx] <- as.numeric(minutes_matches[minutes_idx, "minutes"]) / 60
  }

  result
}
