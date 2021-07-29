#' Adds a parameter for dose intensity
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param tadm_code Total Doses Administered parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total dose administered.  This will be the numerator of the dose
#'   intensity calculation.
#'
#'   Permitted Values: character value
#'
#' @param tpadm_code Total Doses Planned parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the total planned dose.  This will be the denominator of the dose
#'   intensity calculation.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @inheritParams derive_derived_param
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adec adex
#'
#' @export
#'
#' @examples
#' adex <- tibble::tribble(
#` ~USUBJID, ~PARAMCD, ~VISIT, ~ANL01FL, ~ASTDT,            ~AENDT,            ~AVAL, ~AVALU,
#` "P001",   "TNDOSE", "V1",   "Y",      ymd("2020-01-01"), ymd("2020-01-30"), 59,    "mg",
#` "P001",   "TSNDOSE","V1",   "Y",      ymd("2020-01-01"), ymd("2020-02-01"), 96,    "mg",
#` "P001",   "TNDOSE", "V2",   "Y",      ymd("2020-02-01"), ymd("2020-03-15"), 88,    "mg",
#` "P001",   "TSNDOSE","V2",   "Y",      ymd("2020-02-05"), ymd("2020-03-01"), 88,    "mg"
#` )
#`
#` derive_param_doseint(adex,
#`                     filter = NULL,
#`                     new_param = "TNDOSINT",
#`                     tadm_code = "TNDOSE",
#`                     tpadm_code = "TSNDOSE",
#`                     by_vars=vars(USUBJID, VISIT))

derive_param_doseint <- function(dataset,
                                 filter = NULL,
                                 new_param = "TNDOSINT",
                                 tadm_code = "TNDOSE",
                                 tpadm_code = "TSNDOSE",
#                                 drop_values_from=NULL,
                                 by_vars) {

  # TODO Not sure we need to include the drop_values_from if we depend on
  # get_constant_vars in derive_derived_param?
  # Should it be included and set to NULL?

  assert_character_scalar(new_param)
  assert_character_scalar(tadm_code)
  assert_character_scalar(tpadm_code)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset,
                    required_vars = vars(!!!by_vars, PARAMCD, AVAL))

  # TODO Check Units from tadm_code and padm_code are identical if AVALU
  #      exists in dataset, issue warning if they do not match between
  #      parameters

  # TODO Roche specs dictate start (ASTDT) and end (AENDT) dates are taken from
  # another parameter.  Should we allow the user to pass in an additional test
  # code and specify values to retain from that test code or should that be
  # done in another function? Should we leave that to the user to sort out and
  # give examples in the workflow vignette on ways to calculate using
  # Admiral functions?

  derive_derived_param(dataset,
                       filter = !!filter,
                       parameters = c(tadm_code, tpadm_code),
                       by_vars = by_vars,
                       analysis_value = (!!sym(paste0("AVAL.", tadm_code)) / !!sym(paste0("AVAL.", tpadm_code))*100),
                       set_values_to = vars(PARAMCD = !!new_param,
                                            AVALU = "%")
                       #drop_values_from = drop_values_from
                       )
}
