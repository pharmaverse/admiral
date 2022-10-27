admiral_env <- new.env(parent = emptyenv())
admiral_env$admiral_options <- list()

#' Get the value of an admiral option
#'
#' Get the value of an admiral option which can be modified for advanced users.
#'
#' @param option An unquoted expression of commonly used admiral function inputs.
#'
#'   As of now, support only available for `r enumerate(names(admiral_env$admiral_options))`.
#'   See `set_admiral_options()` for a description of the options.
#'
#' @details
#' This function allows flexibility for function inputs that may need to be repeated
#' multiple times in a script, such as `subject_keys`.
#'
#' @author Zelos Zhu
#'
#' @return
#' The value of the specified option unless user changes the default
#' using `set_admiral_options()`.
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
#'   by_vars = get_admiral_option(subject_keys)
#' ) %>%
#'   select(STUDYID, USUBJID, VSTESTCD, VISIT, VSTPT, VSSTRESN, AGE, AGEU)
get_admiral_option <- function(option) {
  # Check for valid option - catch function abuse
  assert_symbol(enquo(option))

  # Find which admiral_options is being called upon
  requested <- as_name(enquo(option))
  possible_inputs <- names(admiral_env$admiral_options)

  if (requested %in% possible_inputs) {
    return(admiral_env$admiral_options[[which(requested == possible_inputs)]])
  }

  # Return message otherwise, catch typos
  else {
    default_err_msg <- sprintf(paste(
      "Invalid function argument, select one unquoted of:",
      enumerate(possible_inputs, quote_fun = dquote, conjunction = "or")
    ))
    abort(default_err_msg)
  }
}

#' Set the values of admiral options
#'
#' Set the values of admiral options that can be modified for advanced users.
#'
#' @param subject_keys Variables to uniquely identify a subject, defaults to
#'   `vars(STUDYID, USUBJID)`. This option is used as default value for the
#'   `subject_keys` argument in all admiral functions.
#'
#' @param set_default If set to `TRUE` will restore admiral defaults.
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
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
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
#' set_admiral_options(set_default = TRUE)
set_admiral_options <- function(subject_keys,
                                # future_input,
                                set_default = FALSE) {
  if (!missing(subject_keys)) {
    assert_vars(subject_keys)
    admiral_env$admiral_options$subject_keys <- subject_keys
  }

  # if (!missing(future_input)) {
  #   assert_vars(future_input) nolint
  #   admiral_env$admiral_options$future_input = future_input nolint
  # }

  assert_logical_scalar(set_default)
  if (set_default == TRUE) {
    admiral_env$admiral_options$subject_keys <- vars(STUDYID, USUBJID)
    # admiral_env$admiral_options$future_input <- vars(...) nolint

    ######################################################
    ### To enhance feature and add inputs as necessary ###
    ######################################################

    # 1. Add additional options such as future_input as shown commented above
    # 2. Update @params with future_input in set_admiral_options roxygen documentation
  }
}
set_admiral_options(set_default = TRUE)
