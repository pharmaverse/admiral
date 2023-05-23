#' Add an Extreme Event Parameter
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_extreme_records()` instead.
#'
#' Add a new parameter for the first or last event occurring in a dataset. The
#' variable given in `new_var` indicates if an event occurred or not. For
#' example, the function can derive a parameter for the first disease
#' progression.
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
#'   `order` argument (if applicable) are expected.
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
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`.
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
#'   A named list returned by `exprs()` defining the variables to be set for the
#'   new parameter, e.g. `exprs(PARAMCD = "PD", PARAM = "Disease Progression")`
#'   is expected. The values must be symbols, character strings, numeric values,
#'   `NA`, or an expression. Note, if you require a date or datetime variable to
#'   be populated, this needs to be defined here.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of symbols created using `exprs()` is expected.
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the
#'   observations of the source dataset (`dataset_source`) restricted by
#'   `filter_source` are not unique with respect to the subject keys
#'   (`subject_key` argument) and `order`.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details
#'   1. The source dataset (`dataset_source`) is restricted to observations fulfilling
#'   `filter_source`.
#'   1. For each subject (with respect to the variables specified for the
#'   `subject_keys` argument) either the first or last observation from the restricted
#'   source dataset is selected. This is depending on `mode`, (with respect to `order`,
#'   if applicable) where the event condition (`filter_source` argument) is fulfilled.
#'   1. For each observation in `dataset_adsl` a new observation is created. For
#'   subjects with event `new_var` is set to `true_value`. For all other
#'   subjects `new_var` is set to `false_value`.
#'   For subjects with event all variables from `dataset_source` are kept. For
#'   subjects without event all variables which are in both `dataset_adsl` and
#'   `dataset_source` are kept.
#'   1. The variables specified by the `set_values_to` argument are added to
#'   the new observations.
#'   1. The new observations are added to input dataset.
#'
#'
#' @return The input dataset with a new parameter indicating if and when an
#'   event occurred
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
derive_param_extreme_event <- function(dataset = NULL,
                                       dataset_adsl,
                                       dataset_source,
                                       filter_source,
                                       order = NULL,
                                       new_var = NULL,
                                       true_value = "Y",
                                       false_value = "N",
                                       mode = "first",
                                       subject_keys = get_admiral_option("subject_keys"),
                                       set_values_to,
                                       check_type = "warning") {
  deprecate_warn("0.11.0", "derive_param_extreme_event()", "derive_extreme_records()")

  # Check input arguments
  filter_source <- assert_filter_cond(enexpr(filter_source))
  assert_vars(subject_keys)
  assert_expr_list(order, optional = TRUE)
  assert_data_frame(dataset_source,
    required_vars = exprs(!!!subject_keys, !!!extract_vars(order))
  )
  new_var <- assert_symbol(enexpr(new_var), optional = TRUE)
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
  if (!is.null(set_values_to$PARAMCD) && !is.null(dataset)) {
    assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  }

  derive_extreme_records(
    dataset,
    dataset_add = dataset_source,
    dataset_ref = dataset_adsl,
    by_vars = subject_keys,
    order = order,
    mode = mode,
    filter_add = !!filter_source,
    check_type = check_type,
    exist_flag = !!new_var,
    true_value = true_value,
    false_value = false_value,
    set_values_to = set_values_to
  )
}
