#' Get the Value of an Admiral Option
#'
#' Get the Value of an Admiral Option Which Can Be Modified for Advanced Users.
#'
#' @param option A character scalar of commonly used admiral function inputs.
#'
#'   As of now, support only available for `r enumerate(names(admiral_environment$admiral_options), quote_fun = dquote, conjunction = "or")`.
#'   See `set_admiral_options()` for a description of the options.
#'
#' @details
#' This function allows flexibility for function inputs that may need to be repeated
#' multiple times in a script, such as `subject_keys`.
#'
#' @author Zelos Zhu
#'
#' @return
#' The value of the specified option.
#'
#' @keywords admiral_options
#' @family admiral_options
#'
#' @export
#'
#' @seealso [vars()], [set_admiral_options()], [derive_param_exist_flag()],
#' [derive_param_first_event()], [derive_param_tte()], [derive_var_disposition_status()],
#' [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_disposition_reason()],
#' [derive_vars_period()], [create_period_dataset()]
#'
#'
#' @examples
#' library(admiral.test)
#' library(dplyr, warn.conflicts = FALSE)
#' data("admiral_vs")
#' data("admiral_dm")
#'
#' # Merging all dm variables to vs
#' derive_vars_merged(
#'   admiral_vs,
#'   dataset_add = select(admiral_dm, -DOMAIN),
#'   by_vars = get_admiral_option("subject_keys")
#' ) %>%
#'   select(STUDYID, USUBJID, VSTESTCD, VISIT, VSTPT, VSSTRESN, AGE, AGEU)
get_admiral_option <- function(option) {
  # Check for valid option - catch function abuse
  assert_character_scalar(option)

  # Find which admiral_options is being called upon
  possible_inputs <- names(admiral_environment$admiral_options)

  if (option %in% possible_inputs) {
    return(admiral_environment$admiral_options[[option]])
  }

  # Return message otherwise, catch typos
  err_msg <- paste(
    "Invalid function argument, select one of:",
    enumerate(possible_inputs, quote_fun = dquote, conjunction = "or")
  )
  abort(err_msg)
}

#' Set the Value of Admiral Options
#'
#' Set the Values of Admiral Options That Can Be Modified for Advanced Users.
#'
#' @param subject_keys Variables to uniquely identify a subject, defaults to
#'   `vars(STUDYID, USUBJID)`. This option is used as default value for the
#'   `subject_keys` argument in all admiral functions.
#'
#' @details
#' Modify an admiral option, e.g `subject_keys`, such that it automatically affects downstream
#' function inputs where `get_admiral_option()` is called such as `derive_param_exist_flag()`.
#'
#' @author Zelos Zhu
#'
#' @return
#' No return value, called for side effects.
#'
#' @keywords admiral_options
#' @family admiral_options
#'
#' @export
#'
#' @seealso [vars()], [get_admiral_option()], [derive_param_exist_flag()],
#' [derive_param_first_event()], [derive_param_tte()], [derive_var_disposition_status()],
#' [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_disposition_reason()],
#' [derive_vars_period()], [create_period_dataset()]
#'
#' @examples
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#' library(tibble)
#' set_admiral_options(subject_keys = vars(STUDYID, USUBJID2))
#'
#' # Derive a new parameter for measurable disease at baseline
#' adsl <- tribble(
#'   ~USUBJID2,
#'   "1",
#'   "2",
#'   "3"
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' tu <- tribble(
#'   ~USUBJID2, ~VISIT,      ~TUSTRESC,
#'   "1",       "SCREENING", "TARGET",
#'   "1",       "WEEK 1",    "TARGET",
#'   "1",       "WEEK 5",    "TARGET",
#'   "1",       "WEEK 9",    "NON-TARGET",
#'   "2",       "SCREENING", "NON-TARGET",
#'   "2",       "SCREENING", "NON-TARGET"
#' ) %>%
#'   mutate(
#'     STUDYID = "XX1234",
#'     TUTESTCD = "TUMIDENT"
#'   )
#'
#' derive_param_exist_flag(
#'   dataset_adsl = adsl,
#'   dataset_add = tu,
#'   filter_add = TUTESTCD == "TUMIDENT" & VISIT == "SCREENING",
#'   condition = TUSTRESC == "TARGET",
#'   false_value = "N",
#'   missing_value = "N",
#'   set_values_to = vars(
#'     PARAMCD = "MDIS",
#'     PARAM = "Measurable Disease at Baseline"
#'   )
#' )
set_admiral_options <- function(subject_keys) {
  if (!missing(subject_keys)) {
    assert_vars(subject_keys)
    admiral_environment$admiral_options$subject_keys <- subject_keys
  }

  # Add future input to function formals above
  # if (!missing(future_input)) {
  #   assert_vars(future_input) nolint
  #   admiral_environment$admiral_options$future_input <- future_input nolint
  # }
}
