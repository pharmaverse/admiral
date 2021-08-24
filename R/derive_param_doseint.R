#' Adds a Parameter for Dose Intensity
#'
#' Adds a record for the dose intensity for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' Total Dose / Planned Dose * 100
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and the `unit_var` parameter,
#'   `PARAMCD`, and `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `tadm_code` and `padm_code`.
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param tadm_code Total Doses Administered parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total dose administered.  The `AVAL` associated with this `PARAMCD`
#'   will be the numerator of the dose intensity calculation.
#'
#'   Permitted Values: character value
#'
#' @param tpadm_code Total Doses Planned parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total planned dose.  The `AVAL` associated with this `PARAMCD`
#'   will be the denominator of the dose intensity calculation.
#'
#'   Permitted Values: character value
#'
#' @param zero_doses Flag indicating logic for handling 0 planned or
#' administered doses for a `by_vars` group
#'
#'   If some instances, the planned dose (`tpadm_code`) or administered dose
#'   (`adm_code`) will be `NA` or 0.  The following options may be applied for cases
#'   where the planned dose (`tpadm_code`) and/or (`tadm_code`) are 0.
#'
#'   Permitted Values:

#'   NULL: Default, No record is returned is either the planned (`tpadm_code`) or
#'   administered (`tadm_code`) are `NA`.  If the planned dose (`tpadm_code`) is
#'   0, `Inf` is returned. If the administered dose (`tadm_code`) is 0 and
#'   the planned dose (`tpadm_code`) is > 0, 0 is returned. A record must exist
#'   for both `tadm_code` and `tpadm_code` for the the dose intensity calculation
#'   to be done.
#'
#'   `Y`: Returns 100 when the planned dose (`tpadm_code`) is 0 and the
#'   administered dose (`tadm_code`) is > 0. Returns 0 when the planned dose
#'   (`tpadm_code`) is 0 and the administered dose (`tadm_code`) is 0.
#'
#' @inheritParams derive_derived_param
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with the new parameter rows added
#'
#' @keywords derivation adec adex
#'
#' @export
#'
#' @examples
#' adex <- tibble::tribble(
#` ~USUBJID, ~PARAMCD, ~VISIT, ~ANL01FL, ~ASTDT,            ~AENDT,            ~AVAL,
#` "P001",   "TNDOSE", "V1",   "Y",      ymd("2020-01-01"), ymd("2020-01-30"), 59,
#` "P001",   "TSNDOSE","V1",   "Y",      ymd("2020-01-01"), ymd("2020-02-01"), 96,
#` "P001",   "TNDOSE", "V2",   "Y",      ymd("2020-02-01"), ymd("2020-03-15"), 88,
#` "P001",   "TSNDOSE","V2",   "Y",      ymd("2020-02-05"), ymd("2020-03-01"), 88
#` )
#`
#` derive_param_doseint(adex,
#'                      by_vars=vars(USUBJID, VISIT),
#'                      set_values_to = vars(PARAMCD = "TNDOSINT",
#'                                           PARAM = "Dose Intensity (%)"),
#`                      tadm_code = "TNDOSE",
#`                      tpadm_code = "TSNDOSE")

derive_param_doseint <- function(dataset,
                                 by_vars,
                                 set_values_to = vars(PARAMCD = "TNDOSINT",
                                                      PARAM = "Dose Intensity (%)"),
                                 tadm_code = "TNDOSE",
                                 tpadm_code = "TSNDOSE",
                                 zero_doses = NULL,
                                 filter = NULL) {

  assert_character_scalar(tadm_code)
  assert_character_scalar(tpadm_code)
  assert_character_scalar(zero_doses, values = c(NULL, "Y"), optional=TRUE)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset,
                    required_vars = vars(!!!by_vars, PARAMCD, AVAL)
                    )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  # Create Dose intensity records
  dataset <-
    derive_derived_param(dataset,
                         filter = !!filter,
                         parameters = c(tadm_code, tpadm_code),
                         by_vars = by_vars,
                         analysis_value = (!!sym(paste0("AVAL.", tadm_code)) / !!sym(paste0("AVAL.", tpadm_code))*100),
                         set_values_to = vars(!!!set_values_to,
                                              tmp_planned_dose_flag = "Y",
                                              tmp_planned_dose = !!sym(paste0("AVAL.", tpadm_code)),
                                              temp_admin_dose = !!sym(paste0("AVAL.", tadm_code)))
                       )


  # handle 0 doses planned if needed
  if (!quo_is_null(enquo(zero_doses))) {
    data_filtered_planned <- dataset %>%
      filter_if(filter) %>%
      filter(AVAL == 0 & PARAMCD == !!tpadm_code) %>%
      rename(tmp_planned_dose = AVAL) %>%
      select(!!!by_vars, tmp_planned_dose)

    data_filtered_admin <- dataset %>%
      filter_if(filter) %>%
      filter(PARAMCD == !!tadm_code) %>%
      rename(tmp_adm_dose = AVAL) %>%
      select(!!!by_vars, tmp_adm_dose)

    dataset <- left_join(dataset, data_filtered_planned, by = vars2chr(by_vars)) %>%
      left_join(data_filtered_admin, by = vars2chr(by_vars)) %>%
      mutate(AVAL = case_when((tmp_planned_dose_flag == "Y") &
                                (tmp_planned_dose == 0) &
                                (tmp_adm_dose > 0) ~ 100,
                              (tmp_planned_dose_flag == "Y") &
                                (tmp_planned_dose == 0) &
                                (tmp_adm_dose == 0) ~ 0,
                              TRUE ~ AVAL))
  }


  dataset <- select(dataset, -starts_with("tmp"))

}
