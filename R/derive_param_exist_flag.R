#' Add an Existence Flag Parameter
#'
#' Add a new parameter indicating that a certain event exists in a dataset.
#' `AVALC` and `AVAL` indicate if an event occurred or not. For example, the
#' function can derive a parameter indicating if there is measureable disease at
#' baseline.
#'
#' @param dataset Input dataset
#'
#'   The variables specified for `subject_keys` and the `PARAMCD` variable are
#'   expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `subject_keys` are expected. For each subject
#'   (as defined by `subject_keys`) from the specified dataset (`dataset_adsl`),
#'   the existence flag is calculated and added as a new observation to the
#'   input datasets (`dataset`)
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `subject_keys` parameter are expected.
#'
#'   This dataset is used to check if an event occurred or not. Any observation
#'   in the dataset fulfilling the event condition (`condition`) is considered
#'   as an event.
#'
#' @param condition Event condition
#'
#'   The condition is evaluated at the additional dataset (`dataset_add`).
#'
#'   For all subjects where it evaluates as `TRUE` at least once `AVALC` is set
#'   to the true value (`true_value`) for the new observations.
#'
#'   For all subjects where it evaluates as `FALSE` or `NA` for all observations
#'   `AVALC` is set to the false value (`false_value`).
#'
#'   For all subjects not present in the additional dataset `AVALC` is set to
#'   the missing value (`missing_value`).
#'
#' @param true_value True value
#'
#'   For all subjects with at least one observations in the additional dataset
#'   (`dataset_add`) fulfilling the event condition (`condition`), `AVALC` is
#'   set to the specified value (`true_value`).
#'
#'   *Default*: `"Y"`
#'
#'   *Permitted Value*: A character scalar
#'
#' @param false_value False value
#'
#'   For all subjects with at least one observations in the additional dataset
#'   (`dataset_add`) but none of them is fulfilling the event condition
#'   (`condition`), `AVALC` is set to the specified value (`false_value`).
#'
#'   *Default*: `NA_character_`
#'
#'   *Permitted Value*: A character scalar
#'
#' @param missing_value Values used for missing information
#'
#'   For all subjects without an observation in the additional dataset
#'   (`dataset_add`), `AVALC` is set to the specified value (`missing_value`).
#'
#'   *Default*: `NA_character_`
#'
#'   *Permitted Value*: A character scalar
#'
#' @param filter_add Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the parameter is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
#'
#' @param aval_fun Function to map character analysis value (`AVALC`) to numeric
#'   analysis value (`AVAL`)
#'
#'   The (first) argument of the function must expect a character vector and the
#'   function must return a numeric vector.
#'
#'   *Default:* `yn_to_numeric` (see `yn_to_numeric()` for details)
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `vars()` defining the variables to be set for the
#'   new parameter, e.g. `vars(PARAMCD = "MDIS", PARAM = "Measurable Disease at
#'   Baseline")` is expected. The values must be symbols, character strings,
#'   numeric values, or `NA`.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `vars()` is expected.
#'
#' @details
#'   1. The additional dataset (`dataset_add`) is restricted to the observations
#'   matching the `filter_add` condition.
#'
#'   1. For each subject in `dataset_adsl` a new observation is created.
#'
#'       - The `AVALC` variable is added and set to the true value (`true_value`)
#'         if for the subject at least one observation exists in the (restricted)
#'         additional dataset where the condition evaluates to `TRUE`.
#'
#'       - It is set to the false value (`false_value`) if for the subject at least
#'         one observation exists and for all observations the condition evaluates
#'         to `FALSE` or `NA`.
#'
#'       - Otherwise, it is set to the missing value (`missing_value`), i.e., for
#'         those subject not in `dataset_add`.
#'
#'   1. The `AVAL` variable is added and set to `aval_fun(AVALC)`.
#'
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the new observations.
#'
#'   1. The new observations are added to input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with a new parameter indicating if an event
#'   occurred (`AVALC`, `AVAL`, and the variables specified by `subject_keys`
#'   and `set_value_to` are populated for the new parameter)
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' # Derive a new parameter for measurable disease at baseline
#' adsl <- tribble(
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3"
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' tu <- tribble(
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
#'   )
#' )
derive_param_exist_flag <- function(dataset = NULL,
                                    dataset_adsl,
                                    dataset_add,
                                    condition,
                                    true_value = "Y",
                                    false_value = NA_character_,
                                    missing_value = NA_character_,
                                    filter_add = NULL,
                                    aval_fun = yn_to_numeric,
                                    subject_keys = get_admiral_option("subject_keys"),
                                    set_values_to) {
  # Check input parameters
  condition <- assert_filter_cond(enquo(condition))
  assert_character_scalar(true_value)
  assert_character_scalar(false_value)
  assert_character_scalar(missing_value)
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_function(aval_fun)
  assert_vars(subject_keys)
  assert_data_frame(
    dataset,
    required_vars = vars(PARAMCD, !!!subject_keys),
    optional = TRUE
  )
  assert_data_frame(dataset_adsl, required_vars = subject_keys)
  assert_data_frame(dataset_add, required_vars = subject_keys)
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  if (!is.null(dataset)) {
    assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  }

  # Create new observations
  new_obs <- derive_var_merged_exist_flag(
    dataset_adsl,
    dataset_add = dataset_add,
    filter_add = !!filter_add,
    condition = !!condition,
    by_vars = subject_keys,
    new_var = AVALC,
    true_value = true_value,
    false_value = false_value,
    missing_value = missing_value
  )
  new_obs <- call_user_fun(mutate(new_obs, AVAL = aval_fun(AVALC)))

  if (!is.numeric(new_obs$AVAL)) {
    abort(paste(
      "Calling `aval_fun(AVALC)` did not result in a numeric vector.\n",
      "A", typeof(new_obs$AVAL), "vector was returned."
    ))
  }

  new_obs <- mutate(new_obs, !!!set_values_to)

  # Create output dataset
  bind_rows(dataset, new_obs)
}
