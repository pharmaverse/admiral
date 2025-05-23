#' Adds a Parameter for Dose Intensity
#'
#' @description Adds a record for the dose intensity for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' **Note:** This is a wrapper function for the more generic `derive_param_computed()`.
#'
#' The analysis value of the new parameter is derived as
#' Total Dose / Planned Dose * 100
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `tadm_code` and `padm_code`.
#'
#' @param by_vars Grouping variables
#'
#'   Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param tadm_code Total Doses Administered parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total dose administered.  The `AVAL` associated with this `PARAMCD`
#'   will be the numerator of the dose intensity calculation.
#'
#' @permitted character value
#'
#' @param tpadm_code Total Doses Planned parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total planned dose.  The `AVAL` associated with this `PARAMCD`
#'   will be the denominator of the dose intensity calculation.
#'
#' @permitted character value
#'
#' @param zero_doses Flag indicating logic for handling 0 planned or
#' administered doses for a `by_vars` group
#'
#'
#' @permitted `Inf`, `100`
#'
#' No record is returned if either the planned (`tpadm_code`) or administered
#' (`tadm_code`) `AVAL` are `NA`.  No record is returned is a record does not
#' exist for both `tadm_code` and `tpadm_code` for the specified `by_var`.
#'
#' If `zero_doses` = `Inf`:
#'   1. If the planned dose (`tpadm_code`) is 0 and administered dose
#'   (`tadm_code`) is 0, `NaN` is returned.
#'   2. If the planned dose (`tpadm_code`) is 0 and the administered dose
#'   (`tadm_code`) is > 0, `Inf` is returned.
#'
#' If `zero_doses` = `100` :
#'   1. If the planned dose (`tpadm_code`) is 0 and administered dose
#'   (`tadm_code`) is 0, 0 is returned.
#'   2. If the planned dose (`tpadm_code`) is 0 and the administered dose
#'   (`tadm_code`) is > 0, 100 is returned.
#'
#' @inheritParams derive_param_map
#'
#' @inheritParams derive_param_computed
#'
#'
#' @return The input dataset with the new parameter rows added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' adex <- tribble(
#'   ~USUBJID, ~PARAMCD, ~VISIT, ~ANL01FL, ~ASTDT, ~AENDT, ~AVAL,
#'   "P001", "TNDOSE", "V1", "Y", ymd("2020-01-01"), ymd("2020-01-30"), 59,
#'   "P001", "TSNDOSE", "V1", "Y", ymd("2020-01-01"), ymd("2020-02-01"), 96,
#'   "P001", "TNDOSE", "V2", "Y", ymd("2020-02-01"), ymd("2020-03-15"), 88,
#'   "P001", "TSNDOSE", "V2", "Y", ymd("2020-02-05"), ymd("2020-03-01"), 88,
#'   "P002", "TNDOSE", "V1", "Y", ymd("2021-01-01"), ymd("2021-01-30"), 0,
#'   "P002", "TSNDOSE", "V1", "Y", ymd("2021-01-01"), ymd("2021-02-01"), 0,
#'   "P002", "TNDOSE", "V2", "Y", ymd("2021-02-01"), ymd("2021-03-15"), 52,
#'   "P002", "TSNDOSE", "V2", "Y", ymd("2021-02-05"), ymd("2021-03-01"), 0
#' )
#'
#' derive_param_doseint(
#'   adex,
#'   by_vars = exprs(USUBJID, VISIT),
#'   set_values_to = exprs(PARAMCD = "TNDOSINT"),
#'   tadm_code = "TNDOSE",
#'   tpadm_code = "TSNDOSE"
#' )
#'
#' derive_param_doseint(
#'   adex,
#'   by_vars = exprs(USUBJID, VISIT),
#'   set_values_to = exprs(PARAMCD = "TDOSINT2"),
#'   tadm_code = "TNDOSE",
#'   tpadm_code = "TSNDOSE",
#'   zero_doses = "100"
#' )
derive_param_doseint <- function(dataset,
                                 by_vars,
                                 set_values_to = exprs(PARAMCD = "TNDOSINT"),
                                 tadm_code = "TNDOSE",
                                 tpadm_code = "TSNDOSE",
                                 zero_doses = "Inf",
                                 filter = NULL) {
  assert_character_scalar(tadm_code)
  assert_character_scalar(tpadm_code)
  assert_character_scalar(zero_doses, values = c("Inf", "100"), optional = TRUE)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_data_frame(dataset,
    required_vars = exprs(!!!by_vars, PARAMCD, AVAL)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)

  # Create Dose intensity records
  aval_tadm <- sym(paste0("AVAL.", tadm_code))
  aval_tpdm <- sym(paste0("AVAL.", tpadm_code))

  analysis_value <- exprs(AVAL = !!aval_tadm / !!aval_tpdm * 100)

  # handle 0 doses planned if needed
  if (zero_doses == "100") {
    update_aval <- exprs(
      AVAL = case_when(
        !!aval_tpdm == 0 &
          !!aval_tadm > 0 ~ 100,
        !!aval_tpdm == 0 &
          !!aval_tadm == 0 ~ 0,
        TRUE ~ AVAL
      )
    )
  } else {
    update_aval <- NULL
  }

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(tadm_code, tpadm_code),
    by_vars = by_vars,
    set_values_to = expr_c(set_values_to, analysis_value, update_aval)
  )
}
