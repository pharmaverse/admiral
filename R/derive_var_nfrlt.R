#' Derive Nominal Relative Time from First Dose (NFRLT)
#'
#' Derives nominal/planned time from first dose in hours by combining visit day
#' information with timepoint descriptions. The function converts timepoint
#' strings to hours using `convert_xxtpt_to_hours()` and adds them to the
#' day-based offset.
#'
#' @param dataset Input dataset containing timepoint and visit day variables.
#'
#' @permitted A data frame or tibble
#'
#' @param new_var Name of the new variable to create (unquoted). Default is
#'   `NFRLT`.
#'
#' @permitted Unquoted variable name
#'
#' @param tpt_var Timepoint variable containing descriptions like "Pre-dose",
#'   "1H Post-dose", etc. (unquoted).
#'
#' @permitted Unquoted variable name
#'
#' @param visit_day Visit day variable (unquoted). This should be the planned/
#'   nominal visit day (e.g., `VISITDY`).
#'
#' @permitted Unquoted variable name
#'
#' @param first_dose_day The day number considered as the first dose day.
#'   Default is 1. For multiple-dose studies, this is typically Day 1.
#'
#' @permitted Numeric scalar (positive integer)
#'
#' @param infusion_duration Numeric value specifying the duration of infusion
#'   in hours. Passed to `convert_xxtpt_to_hours()`. Default is 1 hour.
#'
#' @permitted Numeric scalar (positive)
#'
#' @details
#' The nominal relative time is calculated as:
#'
#' `NFRLT = (visit_day - first_dose_day) * 24 + timepoint_hours`
#'
#' Where:
#' * `visit_day` is the planned visit day (e.g., 1, 2, 8, 15)
#' * `first_dose_day` is the day of first dose (typically 1)
#' * `timepoint_hours` is derived from the timepoint description using
#'   `convert_xxtpt_to_hours()`
#'
#' **Common Use Cases:**
#'
#' * **Single dose study**: Day 1 only, with samples at various timepoints
#'   (e.g., Pre-dose, 1H, 2H, 4H, 8H, 24H)
#' * **Multiple dose study**: Dosing on multiple days (e.g., Day 1, Day 8,
#'   Day 15) with samples around each dose
#' * **Steady state study**: Multiple daily doses with sampling on specific
#'   days
#'
#' **Important Notes:**
#'
#' * The function assumes `visit_day` represents the nominal/planned day, not
#'   the actual study day
#' * For crossover studies, consider deriving NFRLT separately per period
#' * `NA` values in either `visit_day` or `tpt_var` will result in `NA` for
#'   NFRLT
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
#' # Single dose study - Day 1 only
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
#' # Custom infusion duration (2 hours)
#' adpc_inf <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    1,        "Pre-dose",
#'   "001",    1,        "EOI",
#'   "001",    1,        "1H Post EOI"
#' )
#'
#' derive_var_nfrlt(
#'   adpc_inf,
#'   new_var = NFRLT,
#'   tpt_var = PCTPT,
#'   visit_day = VISITDY,
#'   infusion_duration = 2
#' )
#'
#' # Custom first dose day (e.g., dose on Day 3)
#' adpc_custom <- tribble(
#'   ~USUBJID, ~VISITDY, ~PCTPT,
#'   "001",    3,        "Pre-dose",
#'   "001",    3,        "2H Post-dose",
#'   "001",    3,        "24H Post-dose"
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
                             tpt_var,
                             visit_day,
                             first_dose_day = 1,
                             infusion_duration = 1,
                             range_method = "midpoint") {
  new_var <- assert_symbol(enexpr(new_var))
  tpt_var <- assert_symbol(enexpr(tpt_var))
  visit_day <- assert_symbol(enexpr(visit_day))

  assert_data_frame(dataset, required_vars = exprs(!!tpt_var, !!visit_day))
  assert_numeric_vector(first_dose_day, length = 1)
  assert_numeric_vector(infusion_duration, length = 1)

  # Validate first_dose_day is positive
  if (first_dose_day <= 0) {
    cli_abort(
      "{.arg first_dose_day} must be positive, but is {first_dose_day}."
    )
  }

  # Validate infusion_duration is positive
  if (infusion_duration <= 0) {
    cli_abort(
      "{.arg infusion_duration} must be positive, but is {infusion_duration}."
    )
  }

  # Convert timepoint to hours
  tpt_hours <- convert_xxtpt_to_hours(
    dataset[[as_name(tpt_var)]],
    infusion_duration = infusion_duration,
    range_method = range_method
  )

  # Calculate NFRLT
  dataset %>%
    mutate(
      !!new_var := (!!visit_day - first_dose_day) * 24 + tpt_hours
    )
}
