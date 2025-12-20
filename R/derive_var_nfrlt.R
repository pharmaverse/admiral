#' Derive Nominal Relative Time from First Dose (NFRLT)
#'
#' Derives nominal/planned time from first dose in hours by combining visit day
#' information with timepoint descriptions. The function converts timepoint 
#' strings to hours using `convert_xxtpt_to_hours()` and adds them to the 
#' day-based offset.
#'
#' @param dataset Input dataset containing visit day variable and optionally 
#'   timepoint variable.
#'
#' @permitted A data frame or tibble
#'
#' @param new_var Name of the new variable to create (unquoted). Default is 
#'   `NFRLT`.
#'
#' @permitted Unquoted variable name
#'
#' @param tpt_var Timepoint variable containing descriptions like "Pre-dose", 
#'   "1H Post-dose", etc. (unquoted). If not provided or if the variable doesn't 
#'   exist in the dataset, only the visit day offset is calculated (timepoint 
#'   contribution is 0).
#'
#' @permitted Unquoted variable name (optional)
#'
#' @param visit_day Visit day variable (unquoted). This should be the planned/
#'   nominal visit day (e.g., `VISITDY`). Records with `NA` in this variable 
#'   will have NFRLT set to `NA`.
#'
#' @permitted Unquoted variable name
#'
#' @param first_dose_day The day number considered as the first dose day. 
#'   Default is 1. For multiple-dose studies, this is typically Day 1.
#'
#' @permitted Numeric scalar (positive integer)
#'
#' @param treatment_duration Duration of treatment in hours. Can be either:
#'   * A numeric scalar (used for all records), or
#'   * An unquoted variable name from the dataset (e.g., `EXDUR`) where each 
#'     record can have a different treatment duration
#'   
#'   Passed to `convert_xxtpt_to_hours()`. Must be non-negative. Default is 0 
#'   hours (for instantaneous treatments like oral medications).
#'
#' @permitted Numeric scalar or unquoted variable name (non-negative)
#'
#' @param range_method Method for converting time ranges to single values. 
#'   Options are "midpoint" (default), "start", or "end". Passed to 
#'   `convert_xxtpt_to_hours()`. For example, "0-6h" with midpoint returns 3, 
#'   with start returns 0, with end returns 6.
#'
#' @permitted Character scalar ("midpoint", "start", or "end")
#'
#' @param set_values_to_na An optional condition that marks derived NFRLT values 
#'   as `NA`. For example, `set_values_to_na = VISIT == "UNSCHEDULED"` will set 
#'   NFRLT to `NA` for all unscheduled visits. Can use any variables in the 
#'   dataset.
#'
#' @permitted Condition (optional)
#'
#' @details
#' The nominal relative time is calculated as:
#' 
#' `NFRLT = day_offset * 24 + timepoint_hours`
#' 
#' Where:
#' * `day_offset` is calculated from `visit_day` and `first_dose_day`, accounting 
#'   for the absence of Day 0 in clinical trial convention
#' * `timepoint_hours` is derived from the timepoint description using 
#'   `convert_xxtpt_to_hours()`, or 0 if `tpt_var` is not provided
#'
#' **Handling "No Day 0":**
#' 
#' In clinical trials, day numbering typically follows the convention: 
#' ..., Day -2, Day -1, Day 1, Day 2, ... (no Day 0). This function accounts 
#' for this by adjusting the day offset when `visit_day` is negative and 
#' `first_dose_day` is positive.
#' 
#' For example, with `first_dose_day = 1`:
#' * Day -1 → -24 hours (1 day before Day 1, not 2 days)
#' * Day -7 → -168 hours (7 days before Day 1)
#' * Day 1 → 0 hours (first dose day)
#' * Day 8 → 168 hours (7 days after Day 1)
#' 
#' With `first_dose_day = 7`:
#' * Day -1 → -168 hours (7 days before Day 7, accounting for no Day 0)
#' * Day 1 → -144 hours (6 days before Day 7)
#' * Day 6 → -24 hours (1 day before Day 7)
#' * Day 7 → 0 hours (first dose day)
#'
#' **Common Use Cases:**
#' 
#' * **Single dose study**: Day 1 only, with samples at various timepoints 
#'   (e.g., Pre-dose, 1H, 2H, 4H, 8H, 24H)
#' * **Multiple dose study**: Dosing on multiple days (e.g., Day 1, Day 8, 
#'   Day 15) with samples around each dose
#' * **Screening visits**: Negative visit days (e.g., Day -14, Day -7) before 
#'   first dose
#' * **Steady state study**: Multiple daily doses with sampling on specific 
#'   days
#' * **Oral medications**: Use default `treatment_duration = 0` for 
#'   instantaneous absorption
#' * **IV infusions**: Specify `treatment_duration` as infusion duration in 
#'   hours (scalar) or as a variable name containing duration per record
#' * **Exposure records (EX)**: Can be called without `tpt_var` to derive NFRLT 
#'   based only on visit day
#' * **Unscheduled visits**: Use `set_values_to_na` to set NFRLT to `NA` for 
#'   unscheduled or early discontinuation visits
#' * **Variable treatment durations**: Use a variable name (e.g., `EXDUR`) when 
#'   different subjects or visits have different treatment durations
#'
#' **Important Notes:**
#' 
#' * The function assumes `visit_day` represents the nominal/planned day, not 
#'   the actual study day
#' * Day numbering follows clinical trial convention with no Day 0
#' * For timepoints that span multiple days (e.g., "24H Post-dose"), ensure 
#'   `visit_day` is set to the day when the sample was taken. For example, if 
#'   dosing occurs on Day 3, a "24H Post-dose" sample taken on Day 4 should 
#'   have `visit_day = 4`.
#' * For crossover studies, consider deriving NFRLT separately per period
#' * `NA` values in `visit_day` will automatically result in `NA` for NFRLT 
#'   (no need to use `set_values_to_na` for this case)
#' * `NA` values in `tpt_var` will result in `NA` for NFRLT
#' * `NA` values in the `treatment_duration` variable (if using a variable) will 
#'   result in `NA` for NFRLT for those records
#' * Use `set_values_to_na` when you need to set NFRLT to `NA` based on other 
#'   variables (e.g., `VISIT == "UNSCHEDULED"`), especially when `visit_day` 
#'   is populated but should not be used for the NFRLT calculation
#' * If `tpt_var` is not provided or doesn't exist in the dataset, timepoint 
#'   contribution is assumed to be 0 hours
#'
#' **Setting Special Values:**
#' 
#' If you need to set NFRLT to a specific value (e.g., 99999) for certain visits 
#' instead of `NA`, use `set_values_to_na` first to set them to `NA`, then use 
#' a subsequent `mutate()` call to replace those `NA` values:
#' 
#' ```r
#' dataset %>%
#'   derive_var_nfrlt(
#'     ...,
#'     set_values_to_na = VISIT == "UNSCHEDULED"
#'   ) %>%
#'   mutate(NFRLT = if_else(is.na(NFRLT) & VISIT == "UNSCHEDULED", 99999, NFRLT))
#' ```
#'
#' @return The input dataset with the new nominal relative time variable added.
#'
#' @keywords der_bds_findings
#' @family der_bds_findings
#'
#' @seealso [convert_xxtpt_to_hours()]
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' # Single dose study - Day 1 only (oral medication)
#' adpc <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "1H Post-dose",
#'   "001",    1,        "2H Post-dose",
#'   "001",    1,        "4H Post-dose",
#'   "001",    1,        "24H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' # Study with screening visits (negative days)
#' adpc_screen <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    -14,      "Pre-dose",
#'   "001",    -7,       "Pre-dose",
#'   "001",    -1,       "Pre-dose",
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "2H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_screen,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#' # Returns: -336, -168, -24, 0, 2
#' # Note: Day -1 is 24 hours before Day 1 (no Day 0)
#'
#' # Multiple dose study - Days 1, 8, 15
#' adpc_md <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "2H Post-dose",
#'   "001",    8,        "Pre-dose",
#'   "001",    8,        "2H Post-dose",
#'   "001",    15,       "Pre-dose",
#'   "001",    15,       "2H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_md,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' # First dose on Day 7 instead of Day 1
#' adpc_day7 <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    -1,       "Pre-dose",
#'   "001",    1,        "Pre-dose",
#'   "001",    6,        "Pre-dose",
#'   "001",    7,        "Pre-dose",
#'   "001",    8,        "Pre-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_day7,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   first_dose_day = 7
#' )
#' # Returns: -168, -144, -24, 0, 24
#' # Day -1 is 7 days (168 hours) before Day 7
#' # Day 1 is 6 days (144 hours) before Day 7
#'
#' # IV infusion with 2 hour treatment duration (scalar)
#' adpc_inf <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "EOI",
#'   "001",    1,        "1H Post EOI",
#'   "001",    1,        "10MIN PRE EOI"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_inf,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   treatment_duration = 2
#' )
#'
#' # Variable treatment duration - different per subject
#' adpc_var_dur <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,           ~EXDUR,
#'   "001",    1,        "Pre-dose",       1,
#'   "001",    1,        "EOI",            1,
#'   "001",    1,        "1H POST EOI",    1,
#'   "002",    1,        "Pre-dose",       2,
#'   "002",    1,        "EOI",            2,
#'   "002",    1,        "1H POST EOI",    2
#' )
#'
#' derive_var_nfrlt(
#'   adpc_var_dur,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   treatment_duration = EXDUR  # Variable name!
#' )
#'
#' # Exposure records without timepoint variable
#' ex <- tribble(
#'   ~USUBJID, ~VISITDY,
#'   "001",    1,
#'   "001",    8,
#'   "001",    15
#' )
#'
#' derive_var_nfrlt(
#'   ex,
#'   new_var = NFRLT,
#'   visit_day = VISITDY
#' )
#'
#' # With unscheduled visits
#' adpc_unsched <- tribble(
#'   ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
#'   "001",    1,        "VISIT 1",     "Pre-dose",
#'   "001",    1,        "VISIT 1",     "2H Post-dose",
#'   "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
#'   "001",    NA_real_, "UNSCHEDULED", "2H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_unsched,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT == "UNSCHEDULED"
#' )
#'
#' # With early discontinuation
#' adpc_disc <- tribble(
#'   ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
#'   "001",    1,        "VISIT 1",                           "Pre-dose",
#'   "001",    1,        "VISIT 1",                           "2H Post-dose",
#'   "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "Pre-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_disc,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT == "STUDY DRUG EARLY DISCONTINUATION"
#' )
#'
#' # With multiple exclusion criteria
#' adpc_multi <- tribble(
#'   ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
#'   "001",    1,        "VISIT 1",                           "Pre-dose",
#'   "001",    NA_real_, "UNSCHEDULED",                       "Pre-dose",
#'   "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "Pre-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_multi,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT %in% c("UNSCHEDULED", "STUDY DRUG EARLY DISCONTINUATION")
#' )
#'
#' # Setting NFRLT to 99999 for unscheduled visits
#' adpc_unsched_value <- tribble(
#'   ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
#'   "001",    1,        "VISIT 1",     "Pre-dose",
#'   "001",    1,        "VISIT 1",     "2H Post-dose",
#'   "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
#'   "001",    NA_real_, "UNSCHEDULED", "2H Post-dose"
#' )
#'
#' adpc_unsched_value %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLT,
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY,
#'     set_values_to_na = VISIT == "UNSCHEDULED"
#'   ) %>%
#'   mutate(
#'     NFRLT = if_else(is.na(NFRLT) & VISIT == "UNSCHEDULED", 99999, NFRLT)
#'   )
#'
#' # Custom range method (use end of range instead of midpoint)
#' adpc_range <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "0-6h Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_range,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   range_method = "end"
#' )
#'
#' # Using "Before" and "After" terminology
#' adpc_alt <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Before",
#'   "001",    1,        "1H After",
#'   "001",    1,        "2H After"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_alt,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' # Custom first dose day (e.g., dose on Day 3)
#' # Note: 24H Post-dose sample is on Day 4
#' adpc_custom <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    3,        "Pre-dose",
#'   "001",    3,        "2H Post-dose",
#'   "001",    3,        "8H Post-dose",
#'   "001",    4,        "24H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_custom,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   first_dose_day = 3
#' )
derive_var_nfrlt <- function(dataset,
                             new_var = NFRLT,
                             tpt_var = NULL,
                             visit_day,
                             first_dose_day = 1,
                             treatment_duration = 0,
                             range_method = "midpoint",
                             set_values_to_na = NULL) {
  new_var <- assert_symbol(enexpr(new_var))
  tpt_var <- assert_symbol(enexpr(tpt_var), optional = TRUE)
  visit_day <- assert_symbol(enexpr(visit_day))
  set_values_to_na <- assert_filter_cond(enexpr(set_values_to_na), optional = TRUE)
  
  # Check if tpt_var exists in dataset
  has_tpt_var <- !is.null(tpt_var) && as_name(tpt_var) %in% names(dataset)
  
  # Build required_vars based on whether tpt_var exists
  if (has_tpt_var) {
    required_vars <- exprs(!!tpt_var, !!visit_day)
  } else {
    required_vars <- exprs(!!visit_day)
  }
  
  assert_data_frame(dataset, required_vars = required_vars)
  assert_numeric_vector(first_dose_day, length = 1)
  assert_character_vector(range_method, values = c("start", "end", "midpoint"))
  
  # Validate first_dose_day is positive
  if (first_dose_day <= 0) {
    cli_abort(
      "{.arg first_dose_day} must be positive, but is {first_dose_day}."
    )
  }
  
  # Handle treatment_duration - can be variable name or scalar
  treatment_duration_expr <- enexpr(treatment_duration)
  
  if (is.symbol(treatment_duration_expr) && 
      as_name(treatment_duration_expr) %in% names(dataset)) {
    # It's a variable in the dataset
    treatment_duration_vec <- dataset[[as_name(treatment_duration_expr)]]
    assert_numeric_vector(treatment_duration_vec)
  } else {
    # It's a scalar value
    treatment_duration_vec <- eval(treatment_duration_expr)
    assert_numeric_vector(treatment_duration_vec)
  }
  
  # Validate all values are non-negative
  if (any(treatment_duration_vec < 0, na.rm = TRUE)) {
    cli_abort(
      "{.arg treatment_duration} must be non-negative for all values."
    )
  }
  
  # Convert timepoint to hours (or use 0 if no tpt_var)
  if (has_tpt_var) {
    tpt_hours <- convert_xxtpt_to_hours(
      dataset[[as_name(tpt_var)]],
      treatment_duration = treatment_duration_vec,
      range_method = range_method
    )
  } else {
    tpt_hours <- rep(0, nrow(dataset))
  }
  
  # Calculate NFRLT
  # Handle "no Day 0" issue: when both visit_day and first_dose_day are on 
  # opposite sides of 0, add 1 to account for the missing day
  result <- dataset %>%
    mutate(
      day_diff = !!visit_day - first_dose_day,
      # Add 1 if visit_day is negative and first_dose_day is positive
      # (crossing from negative to positive skips Day 0)
      day_diff_adjusted = if_else(
        !!visit_day < 0 & first_dose_day > 0,
        day_diff + 1,
        day_diff
      ),
      !!new_var := day_diff_adjusted * 24 + tpt_hours
    ) %>%
    select(-day_diff, -day_diff_adjusted)
  
  # Set values to NA based on condition
  if (!is.null(set_values_to_na)) {
    result <- result %>%
      mutate(
        !!new_var := if_else(
          !!set_values_to_na,
          NA_real_,
          !!new_var
        )
      )
  }
  
  result
}