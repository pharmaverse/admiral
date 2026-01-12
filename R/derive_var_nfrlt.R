#' Derive Nominal Relative Time from First Dose (NFRLT)
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Derives nominal/planned time from first dose in hours by combining visit day
#' information with timepoint descriptions. The function converts timepoint
#' strings to hours using `convert_xxtpt_to_hours()` and adds them to the
#' day-based offset. Optionally creates a corresponding unit variable.
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
#' @param new_var_unit Name of the unit variable to create (unquoted). If
#'   specified, a character variable will be created containing the unit of
#'   time exactly as provided in `out_unit`. Common CDISC variables are
#'   `FRLTU` (First Dose Relative Time Unit) or `RRLTU` (Reference Relative
#'   Time Unit). If not specified, no unit variable is created.
#'
#' @permitted Unquoted variable name (optional)
#'
#' @param out_unit Unit of time for the output variable. Options are:
#'   * Days: "day", "days", "d"
#'   * Hours: "hour", "hours", "hr", "hrs", "h" (default: "hours")
#'   * Minutes: "minute", "minutes", "min", "mins"
#'   * Weeks: "week", "weeks", "wk", "wks", "w"
#'
#'   Case-insensitive. The internal calculation is performed in hours, then
#'   converted to the specified unit. If `new_var_unit` is specified, it will
#'   contain the value exactly as provided by the user.
#'
#' @permitted Character scalar (see options above)
#'
#' @param tpt_var Timepoint variable containing descriptions like "Pre-dose",
#'   "1H Post-dose", etc. (unquoted). If not provided or if the variable
#'   doesn't exist in the dataset, only the visit day offset is calculated
#'   (timepoint contribution is 0).
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
#' @param set_values_to_na An optional condition that marks derived NFRLT
#'   values as `NA`. For example, `set_values_to_na = VISIT == "UNSCHEDULED"`
#'   will set NFRLT to `NA` for all unscheduled visits. Can use any variables
#'   in the dataset. When `new_var_unit` is specified, the unit variable will
#'   also be set to `NA` for these records.
#'
#' @permitted Condition (optional)
#'
#' @details
#' The nominal relative time is calculated as:
#'
#' `NFRLT = (day_offset * 24 + timepoint_hours) * conversion_factor`
#'
#' Where:
#' * `day_offset` is calculated from `visit_day` and `first_dose_day`,
#'   accounting for the absence of Day 0 in clinical trial convention
#' * `timepoint_hours` is derived from the timepoint description using
#'   `convert_xxtpt_to_hours()`, or 0 if `tpt_var` is not provided
#' * `conversion_factor` is:
#'   - 1 for "hours" (default)
#'   - 1/24 for "days"
#'   - 1/168 for "weeks" (1/24/7)
#'   - 60 for "minutes"
#'
#' If `new_var_unit` is specified, a character variable is created containing
#' the value of `out_unit` exactly as provided by the user. For example:
#' * `out_unit = "hours"` creates unit variable with value "hours"
#' * `out_unit = "HOURS"` creates unit variable with value "HOURS"
#' * `out_unit = "Days"` creates unit variable with value "Days"
#' * `NA` when the corresponding time value is `NA`
#'
#' This matches the behavior of `derive_vars_duration()` and allows
#' consistency when deriving multiple time variables.
#'
#' **Handling "No Day 0":**
#'
#' In clinical trials, day numbering typically follows the convention:
#' ..., Day -2, Day -1, Day 1, Day 2, ... (no Day 0). This function accounts
#' for this by adjusting the day offset when `visit_day` is negative and
#' `first_dose_day` is positive.
#'
#' For example, with `first_dose_day = 1` and different output units:
#' * Day -1, `out_unit = "hours"` -> -24 hours
#' * Day -1, `out_unit = "days"` -> -1 day
#' * Day -1, `out_unit = "weeks"` -> -0.1429 weeks
#' * Day -1, `out_unit = "minutes"` -> -1440 minutes
#' * Day -7 -> -168 hours, -7 days, -1 week, or -10080 minutes
#' * Day 1 -> 0 (in any unit, first dose day)
#' * Day 8 -> 168 hours, 7 days, 1 week, or 10080 minutes
#'
#' With `first_dose_day = 7`:
#' * Day -1 -> -168 hours, -7 days, -1 week, or -10080 minutes
#' * Day 1 -> -144 hours, -6 days, -0.857 weeks, or -8640 minutes
#' * Day 6 -> -24 hours, -1 day, -0.143 weeks, or -1440 minutes
#' * Day 7 -> 0 (in any unit, first dose day)
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
#' * **Exposure records (EX)**: Can be called without `tpt_var` to derive
#'   NFRLT based only on visit day
#' * **Unscheduled visits**: Use `set_values_to_na` to set NFRLT to `NA` for
#'   unscheduled or early discontinuation visits
#' * **Variable treatment durations**: Use a variable name (e.g., `EXDUR`)
#'   when different subjects or visits have different treatment durations
#' * **Hours output**: Use `out_unit = "hours"` (default) for variables like
#'   `NFRLT` with `FRLTU`
#' * **Days output**: Use `out_unit = "days"` for variables like `NFRLTDY`
#'   with `FRLTU`
#' * **Weeks output**: Use `out_unit = "weeks"` for long-term studies with
#'   weekly dosing
#' * **Minutes output**: Use `out_unit = "minutes"` for very short-term PK
#'   studies or when minute precision is needed
#' * **CDISC compliance**: Use `new_var_unit = FRLTU` for first dose relative
#'   time or `new_var_unit = RRLTU` for reference relative time
#' * **Consistency with duration**: Use the same case for `out_unit` across
#'   `derive_vars_duration()` and `derive_var_nfrlt()` to ensure unit
#'   variables match
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
#' * `NA` values in the `treatment_duration` variable (if using a variable)
#'   will result in `NA` for NFRLT for those records
#' * Use `set_values_to_na` when you need to set NFRLT to `NA` based on other
#'   variables (e.g., `VISIT == "UNSCHEDULED"`), especially when `visit_day`
#'   is populated but should not be used for the NFRLT calculation
#' * If `tpt_var` is not provided or doesn't exist in the dataset, timepoint
#'   contribution is assumed to be 0 hours
#' * When using non-hour units, timepoint contributions are still calculated
#'   in hours first (e.g., "2H Post-dose" = 2 hours), then the entire result
#'   is converted to the specified unit
#' * The unit variable (if created) will contain the exact value provided in
#'   `out_unit`, preserving case and format
#'
#' **Setting Special Values:**
#'
#' If you need to set NFRLT to a specific value (e.g., 99999) for certain
#' visits instead of `NA`, use `set_values_to_na` first to set them to `NA`,
#' then use a subsequent `mutate()` call to replace those `NA` values:
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
#' @return The input dataset with the new nominal relative time variable added,
#'   and optionally the unit variable if `new_var_unit` is specified.
#'
#' @keywords der_bds_findings experimental
#' @family der_bds_findings
#'
#' @seealso [convert_xxtpt_to_hours()], [derive_vars_duration()]
#'
#' @export
#'
#' @examplesx
#'
#' @caption Single dose study
#' @info Day 1 only with oral medication
#' @code
#' library(dplyr)
#' library(tibble)
#'
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
#' @caption Single dose study with unit variable
#' @info Creating NFRLT with FRLTU unit variable
#' @code
#' derive_var_nfrlt(
#'   adpc,
#'   new_var = NFRLT,
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Single dose study with different output units
#' @info Deriving NFRLT in different time units with unit variables
#' @code
#' adpc %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLT,
#'     new_var_unit = FRLTU,
#'     out_unit = "HOURS",
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY
#'   ) %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLTDY,
#'     new_var_unit = FRLTDYU,
#'     out_unit = "days",
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY
#'   )
#'
#' @caption Study with screening visits
#' @info Handling negative visit days (no Day 0 in clinical trials)
#' @code
#' adpc_screen <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    -14,      "Screening",
#'   "001",    -7,       "Pre-dose",
#'   "001",    -1,       "Pre-dose",
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "2H Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_screen,
#'   new_var = NFRLT,
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Multiple dose study
#' @info Dosing on Days 1, 8, and 15
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Multiple dose study with days output
#' @info Deriving both NFRLT (hours) and NFRLTDY (days) with unit variables
#' @code
#' adpc_md %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLT,
#'     new_var_unit = FRLTU,
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY
#'   ) %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLTDY,
#'     new_var_unit = FRLTDYU,
#'     out_unit = "days",
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY
#'   )
#'
#' @caption Weekly dosing study
#' @info Long-term study with weekly dosing, using weeks output
#' @code
#' adpc_weekly <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    8,        "Pre-dose",
#'   "001",    15,       "Pre-dose",
#'   "001",    22,       "Pre-dose",
#'   "001",    29,       "Pre-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_weekly,
#'   new_var = NFRLTWK,
#'   new_var_unit = FRLTU,
#'   out_unit = "weeks",
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Short-term PK study with minutes
#' @info Very short timepoints requiring minute precision
#' @code
#' adpc_short <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "5 MIN POST",
#'   "001",    1,        "15 MIN POST",
#'   "001",    1,        "30 MIN POST",
#'   "001",    1,        "1H POST"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_short,
#'   new_var = NFRLTMIN,
#'   new_var_unit = FRLTU,
#'   out_unit = "minutes",
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Custom first dose day
#' @info First dose on Day 7 instead of Day 1
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   first_dose_day = 7
#' )
#'
#' @caption IV infusion with scalar treatment duration
#' @info 2-hour infusion duration for all records
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   treatment_duration = 2
#' )
#'
#' @caption Variable treatment duration
#' @info Different treatment durations per subject using a variable
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   treatment_duration = EXDUR
#' )
#'
#' @caption Exposure records without timepoint variable
#' @info Deriving NFRLT based only on visit day
#' @code
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
#'   new_var_unit = FRLTU,
#'   visit_day = VISITDY
#' )
#'
#' @caption Exposure records with different output units
#' @info Deriving NFRLT in hours, days, and weeks for exposure records
#' @code
#' ex %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLT,
#'     new_var_unit = FRLTU,
#'     visit_day = VISITDY
#'   ) %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLTDY,
#'     new_var_unit = FRLTDYU,
#'     out_unit = "days",
#'     visit_day = VISITDY
#'   ) %>%
#'   derive_var_nfrlt(
#'     new_var = NFRLTWK,
#'     new_var_unit = FRLTWKU,
#'     out_unit = "weeks",
#'     visit_day = VISITDY
#'   )
#'
#' @caption Unscheduled visits
#' @info Setting NFRLT to NA for unscheduled visits
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT == "UNSCHEDULED"
#' )
#'
#' @caption Early discontinuation visits
#' @info Handling study drug early discontinuation
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT == "STUDY DRUG EARLY DISCONTINUATION"
#' )
#'
#' @caption Multiple exclusion criteria
#' @info Excluding multiple visit types
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   set_values_to_na = VISIT %in% c(
#'     "UNSCHEDULED",
#'     "STUDY DRUG EARLY DISCONTINUATION"
#'   )
#' )
#'
#' @caption Setting special values instead of NA
#' @info Using mutate to set NFRLT to 99999 for unscheduled visits
#' @code
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
#'     new_var_unit = FRLTU,
#'     tpt_var = PCTPT,
#'     visit_day = VISITDY,
#'     set_values_to_na = VISIT == "UNSCHEDULED"
#'   ) %>%
#'   mutate(
#'     NFRLT = if_else(is.na(NFRLT) & VISIT == "UNSCHEDULED", 99999, NFRLT),
#'     FRLTU = if_else(is.na(FRLTU) & VISIT == "UNSCHEDULED", "", FRLTU)
#'   )
#'
#' @caption Custom range method
#' @info Using end of range instead of midpoint
#' @code
#' adpc_range <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "0-6h Post-dose"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_range,
#'   new_var = NFRLT,
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   range_method = "end"
#' )
#'
#' @caption Alternative terminology
#' @info Using "Before" and "After" terminology
#' @code
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
#'   new_var_unit = FRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
#'
#' @caption Reference relative time with RRLTU
#' @info Using RRLTU for reference relative time instead of first dose
#' @code
#' derive_var_nfrlt(
#'   adpc,
#'   new_var = NRRLT,
#'   new_var_unit = RRLTU,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   first_dose_day = 8
#' )
#'
#' @caption Case sensitivity in out_unit
#' @info Unit variable preserves the case provided in out_unit
#' @code
#' derive_var_nfrlt(
#'   adpc,
#'   new_var = NFRLT,
#'   new_var_unit = FRLTU,
#'   out_unit = "HOURS",
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY
#' )
derive_var_nfrlt <- function(dataset,
                             new_var = NFRLT,
                             new_var_unit = NULL,
                             out_unit = "HOURS",
                             tpt_var = NULL,
                             visit_day,
                             first_dose_day = 1,
                             treatment_duration = 0,
                             range_method = "midpoint",
                             set_values_to_na = NULL) {
  new_var <- assert_symbol(enexpr(new_var))
  new_var_unit <- assert_symbol(enexpr(new_var_unit), optional = TRUE)
  tpt_var <- assert_symbol(enexpr(tpt_var), optional = TRUE)
  visit_day <- assert_symbol(enexpr(visit_day))
  set_values_to_na <- assert_filter_cond(enexpr(set_values_to_na), optional = TRUE)

  # Store original out_unit before validation (for unit variable)
  original_out_unit <- out_unit

  # Validate out_unit (case-insensitive, returns lowercase)
  out_unit_lower <- assert_character_scalar(
    out_unit,
    values = c(
      c("day", "days", "d"),
      c("hour", "hours", "hr", "hrs", "h"),
      c("minute", "minutes", "min", "mins"),
      c("week", "weeks", "wk", "wks", "w")
    ),
    case_sensitive = FALSE
  )

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

  # Determine conversion factor for unit (using lowercase normalized value)
  conversion_factor <- switch(out_unit_lower,
    day = 1 / 24,
    days = 1 / 24,
    d = 1 / 24,
    hour = 1,
    hours = 1,
    hr = 1,
    hrs = 1,
    h = 1,
    minute = 60,
    minutes = 60,
    min = 60,
    mins = 60,
    week = 1 / 168,
    weeks = 1 / 168,
    wk = 1 / 168,
    wks = 1 / 168,
    w = 1 / 168
  )

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
      !!new_var := (day_diff_adjusted * 24 + tpt_hours) * conversion_factor
    ) %>%
    select(-day_diff, -day_diff_adjusted)

  # Add unit variable if requested (using original case from user)
  if (!is.null(new_var_unit)) {
    result <- result %>%
      mutate(
        !!new_var_unit := if_else(
          is.na(!!new_var),
          NA_character_,
          original_out_unit
        )
      )
  }

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

    # Also set unit variable to NA if it exists
    if (!is.null(new_var_unit)) {
      result <- result %>%
        mutate(
          !!new_var_unit := if_else(
            !!set_values_to_na,
            NA_character_,
            !!new_var_unit
          )
        )
    }
  }

  result
}
