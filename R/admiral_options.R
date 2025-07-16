#' Get the Value of an Admiral Option
#'
#' Get the Value of an Admiral Option Which Can Be Modified for Advanced Users.
#'
#' @param option A character scalar of commonly used admiral function inputs.
#'
#'   As of now, support only available for
#'   `r ansi_collapse(map_chr(names(admiral_environment$admiral_options), ~ paste0('"', ., '"')))`.
#'   See `set_admiral_options()` for a description of the options.
#'
#' @details
#' This function allows flexibility for function inputs that may need to be repeated
#' multiple times in a script, such as `subject_keys`.
#'
#'
#' @return
#' The value of the specified option.
#'
#' @keywords admiral_options
#' @family admiral_options
#'
#' @export
#'
#' @seealso [set_admiral_options()], [derive_param_exist_flag()], [derive_param_tte()]
#'  [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_period()],
#'  [create_period_dataset()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' dm <- tribble(
#'   ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",   "DM", "01-1302",   61, "YEARS",
#'   "PILOT01",   "DM", "17-1344",   64, "YEARS"
#' )
#'
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~VSTESTCD,     ~VISIT,     ~VSTPT, ~VSSTRESN,
#'   "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE",    "LYING",        76,
#'   "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE", "STANDING",        87,
#'   "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2",    "LYING",        71,
#'   "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2", "STANDING",        79,
#'   "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE",    "LYING",        88,
#'   "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE", "STANDING",        86,
#'   "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2",    "LYING",        84,
#'   "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2", "STANDING",        82
#' )
#'
#' # Merging all dm variables to vs
#' derive_vars_merged(
#'   vs,
#'   dataset_add = select(dm, -DOMAIN),
#'   by_vars = get_admiral_option("subject_keys")
#' )
get_admiral_option <- function(option) {
  # Check for valid option - catch function abuse
  assert_character_scalar(option)

  # Find which admiral_options is being called upon
  possible_inputs <- names(admiral_environment$admiral_options)

  if (option %in% possible_inputs) {
    return(admiral_environment$admiral_options[[option]])
  }

  # change cli `.val` to end with OR instead of AND
  divid <- cli_div(theme = list(.val = list("vec-last" = ", or ", "vec_sep2" = " or ")))
  # Return message otherwise, catch typos
  cli_abort(c(
    "Invalid function argument.",
    "i" = "Select one of {.val {possible_inputs}}"
  ))
}

#' Set the Value of admiral Options
#'
#' Set the values of admiral options that can be modified for advanced users.
#'
#' @param subject_keys Variables to uniquely identify a subject, defaults to
#'   `exprs(STUDYID, USUBJID)`. This option is used as default value for the
#'   `subject_keys` argument in all admiral functions.
#'
#' @param signif_digits Holds number of significant digits when comparing to numeric variables,
#'  defaults to `15`. This option is used as default value for the  `signif_dig` argument in
#'  admiral functions `derive_var_atoxgr_dir()` and `derive_var_anrind()`.
#'
#' @param save_memory If set to `TRUE`, an alternative algorithm is used in the
#'   functions `derive_vars_joined()`, `derive_var_joined_exist_flag()`,
#'   `derive_extreme_event()`, and `filter_joined()` which requires less memory
#'   but more run-time.
#'
#' @details
#' Modify an admiral option, e.g `subject_keys`, such that it automatically affects downstream
#' function inputs where `get_admiral_option()` is called such as `derive_param_exist_flag()`.
#'
#'
#' @return
#' No return value, called for side effects.
#'
#' @keywords admiral_options
#' @family admiral_options
#'
#' @export
#'
#' @seealso [get_admiral_option()], [derive_param_exist_flag()],[derive_param_tte()],
#' [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_period()],
#' [create_period_dataset()], [derive_var_atoxgr_dir()], [derive_var_anrind()]
#'
#' @examples
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#' library(tibble)
#' set_admiral_options(subject_keys = exprs(STUDYID, USUBJID2))
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
#'   ~USUBJID2,      ~VISIT,    ~TUSTRESC,
#'   "1",       "SCREENING",     "TARGET",
#'   "1",          "WEEK 1",     "TARGET",
#'   "1",          "WEEK 5",     "TARGET",
#'   "1",          "WEEK 9", "NON-TARGET",
#'   "2",       "SCREENING", "NON-TARGET",
#'   "2",       "SCREENING", "NON-TARGET"
#' ) %>%
#'   mutate(
#'     STUDYID = "XX1234",
#'     TUTESTCD = "TUMIDENT"
#'   )
#'
#' derive_param_exist_flag(
#'   dataset_ref = adsl,
#'   dataset_add = tu,
#'   filter_add = TUTESTCD == "TUMIDENT" & VISIT == "SCREENING",
#'   condition = TUSTRESC == "TARGET",
#'   false_value = "N",
#'   missing_value = "N",
#'   set_values_to = exprs(
#'     PARAMCD = "MDIS",
#'     PARAM = "Measurable Disease at Baseline"
#'   )
#' )
#'
#' set_admiral_options(signif_digits = 14)
#'
#' # Derive ANRIND for ADVS
#' advs <- tribble(
#'   ~PARAMCD, ~AVAL, ~ANRLO, ~ANRHI,
#'   "DIABP",     59,     60,     80,
#'   "SYSBP",    120,     90,    130,
#'   "RESP",      21,      8,     20,
#' )
#'
#' derive_var_anrind(advs)
#'
set_admiral_options <- function(subject_keys, signif_digits, save_memory) {
  if (!missing(subject_keys)) {
    assert_vars(subject_keys)
    admiral_environment$admiral_options$subject_keys <- subject_keys
  }
  if (!missing(signif_digits)) {
    assert_integer_scalar(signif_digits, subset = "positive")
    admiral_environment$admiral_options$signif_digits <- signif_digits
  }
  if (!missing(save_memory)) {
    assert_logical_scalar(save_memory)
    admiral_environment$admiral_options$save_memory <- save_memory
  }

  # Add future input to function formals above
  # if (!missing(future_input)) {
  #   assert_vars(future_input) nolint
  #   admiral_environment$admiral_options$future_input <- future_input nolint
  # }
}
