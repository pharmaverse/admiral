#' Add an Existence Flag Parameter
#'
#' Add a new parameter indicating that a certain event exists in a dataset.
#' `AVALC` and `AVAL` indicate if an event occurred or not. For example, the
#' function can derive a parameter indicating if there is measurable disease at
#' baseline.
#'
#' @param dataset Input dataset
#'
#'   The variables specified for `by_vars` and the `PARAMCD` variable are
#'   expected.
#'
#' @param dataset_ref Reference dataset, e.g., ADSL
#'
#'   The variables specified in `by_vars` are expected. For each group
#'   (as defined by `by_vars`) from the specified dataset (`dataset_ref`),
#'   the existence flag is calculated and added as a new observation to the
#'   input datasets (`dataset`).
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   This dataset is used to check if an event occurred or not. Any observation
#'   in the dataset fulfilling the event condition (`condition`) is considered
#'   as an event.
#'
#' @param condition Event condition
#'
#'   The condition is evaluated at the additional dataset (`dataset_add`).
#'
#'   For all groups where it evaluates as `TRUE` at least once `AVALC` is set
#'   to the true value (`true_value`) for the new observations.
#'
#'   For all groups where it evaluates as `FALSE` or `NA` for all observations
#'   `AVALC` is set to the false value (`false_value`).
#'
#'   For all groups not present in the additional dataset `AVALC` is set to
#'   the missing value (`missing_value`).
#'
#' @param true_value True value
#'
#'   For all groups with at least one observations in the additional dataset
#'   (`dataset_add`) fulfilling the event condition (`condition`), `AVALC` is
#'   set to the specified value (`true_value`).
#'
#'   *Default*: `"Y"`
#'
#'   *Permitted Value*: A character scalar
#'
#' @param false_value False value
#'
#'   For all groups with at least one observations in the additional dataset
#'   (`dataset_add`) but none of them is fulfilling the event condition
#'   (`condition`), `AVALC` is set to the specified value (`false_value`).
#'
#'   *Default*: `NA_character_`
#'
#'   *Permitted Value*: A character scalar
#'
#' @param missing_value Values used for missing information
#'
#'   For all groups without an observation in the additional dataset
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
#'   *Deprecated*, please use `set_values_to` instead.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `exprs()` defining the variables to be set for the
#'   new parameter, e.g. `exprs(PARAMCD = "MDIS", PARAM = "Measurable Disease at
#'   Baseline")` is expected. The values must be symbols, character strings,
#'   numeric values, `NA`, or expressions.
#'
#' @param by_vars Variables to uniquely identify a group
#'
#'   A list of symbols created using `exprs()` is expected.
#'
#' @param dataset_adsl *Deprecated*, please use `dataset_ref` instead.
#'
#' @param subject_keys *Deprecated*, please use `by_vars` instead.
#'
#' @details
#'   1. The additional dataset (`dataset_add`) is restricted to the observations
#'   matching the `filter_add` condition.
#'
#'   1. For each group in `dataset_ref` a new observation is created.
#'
#'       - The `AVALC` variable is added and set to the true value (`true_value`)
#'         if for the group at least one observation exists in the (restricted)
#'         additional dataset where the condition evaluates to `TRUE`.
#'
#'       - It is set to the false value (`false_value`) if for the group at least
#'         one observation exists and for all observations the condition evaluates
#'         to `FALSE` or `NA`.
#'
#'       - Otherwise, it is set to the missing value (`missing_value`), i.e., for
#'         those groups not in `dataset_add`.
#'
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the new observations.
#'
#'   1. The new observations are added to input dataset.
#'
#'
#' @return The input dataset with a new parameter indicating if an event
#'   occurred (`AVALC` and the variables specified by `by_vars`
#'   and `set_value_to` are populated for the new parameter).
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
#'   dataset_ref = adsl,
#'   dataset_add = tu,
#'   filter_add = TUTESTCD == "TUMIDENT" & VISIT == "SCREENING",
#'   condition = TUSTRESC == "TARGET",
#'   false_value = "N",
#'   missing_value = "N",
#'   set_values_to = exprs(
#'     AVAL = yn_to_numeric(AVALC),
#'     PARAMCD = "MDIS",
#'     PARAM = "Measurable Disease at Baseline"
#'   )
#' )
derive_param_exist_flag <- function(dataset = NULL,
                                    dataset_ref,
                                    dataset_add,
                                    condition,
                                    true_value = "Y",
                                    false_value = NA_character_,
                                    missing_value = NA_character_,
                                    filter_add = NULL,
                                    aval_fun,
                                    by_vars = get_admiral_option("subject_keys"),
                                    set_values_to,
                                    dataset_adsl,
                                    subject_keys) {
  ### BEGIN DEPRECATION
  if (!missing(dataset_adsl)) {
    deprecate_warn(
      "0.11.0", "derive_param_exist_flag(dataset_adsl = )",
      "derive_param_exit_flag(dataset_ref = )"
    )
    # assign deprecated argument to new variable
    dataset_ref <- dataset_adsl
  }

  if (!missing(subject_keys)) {
    deprecate_warn(
      "0.11.0", "derive_param_exist_flag(subject_keys = )",
      "derive_param_exit_flag(by_vars = )"
    )
    # assign deprecated argument to new variable
    by_vars <- subject_keys
  }
  ### END DEPRECATION

  # Check input parameters
  condition <- assert_filter_cond(enexpr(condition))
  assert_character_scalar(true_value)
  assert_character_scalar(false_value)
  assert_character_scalar(missing_value)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = exprs(PARAMCD, !!!by_vars),
    optional = TRUE
  )
  assert_data_frame(dataset_ref, required_vars = by_vars)
  assert_data_frame(dataset_add, required_vars = by_vars)
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  if (!is.null(dataset)) {
    assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  }

  if (!missing(aval_fun)) {
    assert_function(aval_fun)
    deprecate_warn(
      "0.11.0",
      "derive_param_exist_flag(aval_fun = )",
      "derive_param_exist_flag(set_values_to = )"
    )
    set_values_to <- exprs(!!!set_values_to, AVAL = aval_fun(AVALC))
  }

  # Create new observations
  new_obs <- derive_var_merged_exist_flag(
    dataset_ref,
    dataset_add = dataset_add,
    filter_add = !!filter_add,
    condition = !!condition,
    by_vars = by_vars,
    new_var = AVALC,
    true_value = true_value,
    false_value = false_value,
    missing_value = missing_value
  )

  new_obs <- process_set_values_to(
    new_obs,
    set_values_to = set_values_to
  )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
