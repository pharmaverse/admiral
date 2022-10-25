admiral_default <- new.env(parent = emptyenv())
admiral_default$admiral_options <- list(
  subject_keys = vars(STUDYID, USUBJID)
  # future_input = vars(...)

  ######################################################
  ### To enhance feature and add inputs as necessary ###
  ######################################################

  # 1. Add additional options such as future_input as shown commented above
  # 2. Update @params with future_input in set_admiral_options roxygen documentation
  # 3. Modify the commented sections of set_admiral_options()
)

#' Call a package-standard input
#'
#' Call a package-standard input that can be modified for advanced users.
#'
#' @param input An unquoted expression of commonly used admiral function inputs.
#'
#'   As of now, support only available for `subject keys`.
#'
#' @details
#' This function allows flexibility for function inputs that may need to be repeated
#' multiple times in a script, such as for `subject keys`.
#'
#' @author Zelos Zhu
#'
#' @return
#' Currently `subject_keys` returns `vars(STUDYID, USUBJID)` unless user changes the default
#' using `set_admiral_options()`.
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @export
#'
#' @seealso [vars()], [set_admiral_options()], [derive_param_exist_flag()],
#' [derive_param_first_event()], [derive_param_tte()], [derive_var_disposition_status()],
#' [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_disposition_reason()]
#'
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' # Derive a new parameter for measurable disease at baseline
#' adsl <- tibble::tribble(
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3"
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' tu <- tibble::tribble(
#'   ~USUBJID, ~VISIT,      ~TUSTRESC,
#'   "1",      "SCREENING", "TARGET",
#'   "1",      "WEEK 1",    "TARGET",
#'   "1",      "WEEK 5",    "TARGET",
#'   "1",      "WEEK 9",    "NON-TARGET",
#'   "2",      "SCREENING", "NON-TARGET",
#'   "2",      "SCREENING", "NON-TARGET"
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
#'   ),
#'   subject_keys = get_admiral_options(subject_keys)
#' )
get_admiral_options <- function(input) {
  # Check for valid input - catch function abuse
  assert_expr(enquo(input))

  # Find which admiral_options is being called upon
  requested <- as_name(enquo(input))
  possible_inputs <- names(admiral_default$admiral_options)

  if (requested %in% possible_inputs) {
    return(admiral_default$admiral_options[[which(requested == possible_inputs)]])
  }

  # Return message otherwise, catch typos
  else {
    default_err_msg <- sprintf(paste(
      "Invalid function argument, select one unquoted of:",
      paste0(possible_inputs, collapse = " or ")
    ))
    abort(default_err_msg)
  }
}

#' Set a package-standard input
#'
#' Set a package-standard input that can be modified for advanced users.
#'
#' @param subject_keys Variables to uniquely identify a subject, defaults to
#'   vars(STUDYID, USUBJID)
#'
#' @param set_default If set to `TRUE` will restore admiral defaults.
#'
#' @details
#' This function allows flexibility for function inputs that may need to be repeated
#' multiple times in a script, such as for `subject keys`.
#'
#' @author Zelos Zhu
#'
#' @return
#' Modify an admiral option, e.g `subject_keys`, such that it automatically affects downstream
#' function inputs where `get_admiral_option()` is called such as `derive_param_exist_flag()`.
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @export
#'
#' @seealso [vars()], [get_admiral_options()], [derive_param_exist_flag()],
#' [derive_param_first_event()], [derive_param_tte()], [derive_var_disposition_status()],
#' [derive_var_dthcaus()], [derive_var_extreme_dtm()], [derive_vars_disposition_reason()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#' set_admiral_options(subject_keys = vars(STUDYID, USUBJID2))
#'
#' # Derive a new parameter for measurable disease at baseline
#' adsl <- tibble::tribble(
#'   ~USUBJID2,
#'   "1",
#'   "2",
#'   "3"
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' tu <- tibble::tribble(
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
    admiral_default$admiral_options$subject_keys <- subject_keys
  }

  # if (!missing(future_input)) {
  #   assert_vars(future_input)
  #   admiral_default$admiral_options$future_input = future_input
  # }

  if (!missing(set_default)) {
    assert_logical_scalar(set_default)
    if (set_default == TRUE) {
      admiral_default$admiral_options$subject_keys <- vars(STUDYID, USUBJID)
      # admiral_default$admiral_options$future_input = vars(...)
    }
  }
}
