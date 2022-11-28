#' Add a First Event Parameter
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_param_extreme_event()` instead with the `order` argument instead of the `date_var` argument.
#'
#' @param dataset Input dataset
#'
#'   The `PARAMCD` variable is expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `subject_keys` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param dataset_source Source dataset
#'
#'   All observations in the specified dataset fulfilling the condition
#'   specified by `filter_source` are considered as event.
#'
#'   The variables specified by the `subject_keys` and
#'   `date_var` parameter are expected.
#'
#' @param filter_source Source filter
#'
#'   All observations in `dataset_source` fulfilling the specified condition are
#'   considered as event.
#'
#'   For subjects with at least one event `AVALC` is set to `"Y"`, `AVAL` to
#'   `1`, and `ADT` to the first date where the condition is fulfilled.
#'
#'   For all other subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to
#'   `NA`.
#'
#' @param date_var Date variable
#'
#'   Date variable in the source dataset (`dataset_source`). The variable is
#'   used to sort the source dataset. `ADT` is set to the specified variable for
#'   events.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `vars()` defining the variables to be set for the
#'   new parameter, e.g. `vars(PARAMCD = "PD", PARAM = "Disease Progression")`
#'   is expected. The values must be symbols, character strings, numeric values,
#'   or `NA`.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `vars()` is expected.
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the
#'   observations of the input dataset restricted to the source parameter
#'   (`source_param`) are not unique with respect to the subject keys
#'   (`subject_key` parameter) and `ADT`.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details
#'   1. The input dataset is restricted to observations fulfilling
#'   `filter_source`.
#'   1. For each subject (with respect to the variables specified for the
#'   `subject_keys` parameter) the first observation (with respect to
#'   `date_var`) where the event condition (`filter_source` parameter) is
#'   fulfilled is selected.
#'   1. For each observation in `dataset_adsl` a new observation is created. For
#'   subjects with event `AVALC` is set to `"Y"`, `AVAL` to `1`, and `ADT` to
#'   the first date where the event condition is fulfilled. For all other
#'   subjects `AVALC` is set to `"N"`, `AVAL` to `0`, and `ADT` to `NA`.
#'   For subjects with event all variables from `dataset_source` are kept. For
#'   subjects without event all variables which are in both `dataset_adsl` and
#'   `dataset_source` are kept.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the new observations.
#'   1. The new observations are added to input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with a new parameter indicating if and when an
#'   event occurred
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
derive_param_first_event <- function(dataset,
                                     dataset_adsl,
                                     dataset_source,
                                     filter_source,
                                     date_var,
                                     subject_keys = vars(STUDYID, USUBJID),
                                     set_values_to,
                                     check_type = "warning") {
  ### DEPRECATION
  deprecate_warn("0.9.0",
    "derive_param_first_event()",
    details = "Please use `derive_param_extreme_event()` instead with the `order` argument instead of the `date_var` argument"
  )

  filter_source <- enquo(filter_source)
  date_var <- enquo(date_var)
  tmp_var <- get_new_tmp_var(dataset = dataset)
  tmp_var <- enquo(tmp_var)

  derive_param_extreme_event(
    dataset = dataset,
    dataset_adsl = dataset_adsl,
    dataset_source = dataset_source,
    filter_source = !!filter_source,
    order = vars(!!date_var),
    new_var = !!tmp_var,
    subject_keys = subject_keys,
    set_values_to = set_values_to,
    check_type = check_type,
    mode = "first"
  ) %>%
    mutate(
      AVALC = coalesce(!!tmp_var, AVALC),
      AVAL = if_else(!!tmp_var == "Y", true = 1, false = 0)
    ) %>%
    remove_tmp_vars()
}

#' Add an Extreme Event Parameter
#'
#' Add a new parameter for the first or last event occurring in a dataset. The
#'  variable given in `new_var` indicates if an event occurred or not. For example,
#'  the function can derive a parameter for the first disease progression.
#'
#' @param dataset Input dataset
#'
#'   The `PARAMCD` variable is expected.
#'
#' @param dataset_adsl ADSL input dataset
#'
#'   The variables specified for `subject_keys` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param dataset_source Source dataset
#'
#'   All observations in the specified dataset fulfilling the condition
#'   specified by `filter_source` are considered as an event.
#'
#'   The variables specified by the `subject_keys` and
#'   `order` parameter (if applicable) are expected.
#'
#' @param filter_source Source filter
#'
#'   All observations in `dataset_source` fulfilling the specified condition are
#'   considered as an event.
#'
#'   For subjects with at least one event `new_var` is set to `true_value`.
#'
#'   For all other subjects `new_var` is set to `false_value`.
#'
#' @param order Order variable
#'
#'   List of symbols for sorting the source dataset (`dataset_source`).
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'   created by `vars()`, e.g., `vars(ADT, desc(AVAL))`.
#'
#' @param new_var New variable
#'
#'   The name of the variable which will indicate whether an event happened or not.
#'
#' @param true_value True value
#'
#'   For all subjects with at least one observation in the source dataset
#'   (`dataset_source`) fulfilling the event condition (`filter_source`),
#'   `new_var` is set to the specified value `true_value`.
#'
#' @param false_value False value
#'
#'   For all other subjects in `dataset_adsl` without an event, `new_var` is set to
#'   the specified value `false_value`.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each subject is selected.
#'   If `"last"` is specified, the last observation of each subject is selected.
#'
#'   *Permitted Values*: `"first"`, `"last"`
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `vars()` defining the variables to be set for the
#'   new parameter, e.g. `vars(PARAMCD = "PD", PARAM = "Disease Progression")`
#'   is expected. The values must be symbols, character strings, numeric values,
#'   or `NA`. Note, if you require a date or datetime variable to be populated,
#'   this needs to be defined here.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `vars()` is expected.
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the
#'   observations of the input dataset restricted to the source parameter
#'   (`source_param`) are not unique with respect to the subject keys
#'   (`subject_key` parameter) and order variables (`order` parameter).
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details
#'   1. The source dataset (`dataset_source`) is restricted to observations fulfilling
#'   `filter_source`.
#'   1. For each subject (with respect to the variables specified for the
#'   `subject_keys` parameter) either the first or last observation from the restricted
#'   source dataset is selected. This is depending on `mode`, (with respect to `order`,
#'   if applicable) where the event condition (`filter_source` parameter) is fulfilled.
#'   1. For each observation in `dataset_adsl` a new observation is created. For
#'   subjects with event `new_var` is set to `true_var`. For all other
#'   subjects `new_var` is set to `false_var`.
#'   For subjects with event all variables from `dataset_source` are kept. For
#'   subjects without event all variables which are in both `dataset_adsl` and
#'   `dataset_source` are kept.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the new observations.
#'   1. The new observations are added to input dataset.
#'
#' @author Stefan Bundfuss Sophie Shapcott
#'
#' @return The input dataset with a new parameter indicating if and when an
#'   event occurred
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' # Derive a new parameter for the first disease progression (PD)
#' adsl <- tribble(
#'   ~USUBJID, ~DTHDT,
#'   "1",      ymd("2022-05-13"),
#'   "2",      ymd(""),
#'   "3",      ymd("")
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' adrs <- tribble(
#'   ~USUBJID, ~ADTC,        ~AVALC,
#'   "1",      "2020-01-02", "PR",
#'   "1",      "2020-02-01", "CR",
#'   "1",      "2020-03-01", "CR",
#'   "1",      "2020-04-01", "SD",
#'   "2",      "2021-06-15", "SD",
#'   "2",      "2021-07-16", "PD",
#'   "2",      "2021-09-14", "PD"
#' ) %>%
#'   mutate(
#'     STUDYID = "XX1234",
#'     ADT = ymd(ADTC),
#'     PARAMCD = "OVR",
#'     PARAM = "Overall Response",
#'     ANL01FL = "Y"
#'   ) %>%
#'   select(-ADTC)
#'
#' derive_param_extreme_event(
#'   adrs,
#'   dataset_adsl = adsl,
#'   dataset_source = adrs,
#'   filter_source = PARAMCD == "OVR" & AVALC == "PD",
#'   order = vars(ADT),
#'   new_var = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   mode = "first",
#'   set_values_to = vars(
#'     PARAMCD = "PD",
#'     PARAM = "Disease Progression",
#'     ANL01FL = "Y",
#'     ADT = ADT
#'   )
#' )
#'
#' # derive parameter indicating death
#' derive_param_extreme_event(
#'   dataset_adsl = adsl,
#'   dataset_source = adsl,
#'   filter_source = !is.na(DTHDT),
#'   new_var = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   mode = "first",
#'   set_values_to = vars(
#'     PARAMCD = "DEATH",
#'     PARAM = "Death",
#'     ANL01FL = "Y",
#'     ADT = DTHDT
#'   )
#' )
derive_param_extreme_event <- function(dataset = NULL,
                                       dataset_adsl,
                                       dataset_source,
                                       filter_source,
                                       order = NULL,
                                       new_var = AVALC,
                                       true_value = "Y",
                                       false_value = "N",
                                       mode = "first",
                                       subject_keys = vars(STUDYID, USUBJID),
                                       set_values_to,
                                       check_type = "warning") {
  # Check input parameters
  filter_source <- assert_filter_cond(enquo(filter_source))
  assert_vars(subject_keys)
  assert_vars(order, optional = TRUE)
  assert_data_frame(dataset_source,
    required_vars = vars(!!!subject_keys, !!!extract_vars(order))
  )
  new_var <- assert_symbol(enquo(new_var))
  assert_same_type(true_value, false_value)
  assert_data_frame(dataset, optional = TRUE)
  assert_data_frame(dataset_adsl, required_vars = subject_keys)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  mode <- assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  if (!is.null(set_values_to$PARAMCD) & !is.null(dataset)) {
    assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  }

  # Create new observations
  source_vars <- colnames(dataset_source)
  adsl_vars <- colnames(dataset_adsl)

  events <- dataset_source %>%
    filter_if(filter_source) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = order,
      mode = mode,
      check_type = check_type
    ) %>%
    mutate(!!new_var := true_value)

  noevents <- anti_join(
    select(dataset_adsl, intersect(source_vars, adsl_vars)),
    select(events, !!!subject_keys),
    by = sapply(subject_keys, as_name)
  ) %>%
    mutate(!!new_var := false_value)

  new_obs <- bind_rows(events, noevents) %>%
    mutate(
      !!!set_values_to
    )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
