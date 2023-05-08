#' Add the First or Last Observation for Each By Group as New Records
#'
#' Add the first or last observation for each by group as new observations. The
#' new observations can be selected from the input dataset or an additional
#' dataset. This function can be used for adding the maximum or minimum value
#' as a separate visit. All variables of the selected observation are kept. This
#' distinguishes `derive_extreme_records()` from `derive_summary_records()`,
#' where only the by variables are populated for the new records.
#'
#' @param dataset Input dataset
#'
#'   If `dataset_add` is not specified, the new records are selected from the
#'   input dataset. In this case the variables specified by `by_vars` and
#'   `order` are expected.
#'
#' @param dataset_ref Reference dataset
#'
#'   The variables specified for `by_vars` are expected. For each
#'   observation of the specified dataset a new observation is added to the
#'   input dataset.
#'
#' @param dataset_add Additional dataset
#'
#'   Observations from the specified dataset are added as new records to the
#'   input dataset (`dataset`).
#'
#'   All observations in the specified dataset fulfilling the condition
#'   specified by `filter_source` are considered. If `mode` and `order` are
#'   specified, the first or last observation within each by group, defined by
#'   `by_vars`, is selected.
#'
#'   If the argument is not specified, the input dataset (`dataset`) is used.
#'
#'   The variables specified by the `by_vars` and `order` argument (if
#'   applicable) are expected.
#'
#' @param by_vars Grouping variables
#'
#'   If `dataset_ref` is specified, this argument must be specified.
#'
#'   *Permitted Values*: list of variables created by `exprs()`
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations in `dataset_add` fulfilling the specified condition are
#'   considered.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is added
#'   to the input dataset. If `"last"` is specified, the last observation of
#'   each by group is added to the input dataset.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) additional dataset are not unique
#'   with respect to the by variables and the order.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @param exist_flag Existence flag
#'
#'   The specified variable is added to the output dataset.
#'
#'   For by groups with at least one observation in the additional dataset
#'   (`dataset_add`) `exist_flag` is set to the value specified by the
#'   `true_value` argument.
#'
#'   For all other by groups `exist_flag` is set to the value specified by the
#'   `false_value` argument.
#'
#'   *Permitted Values:* Variable name
#'
#' @param true_value True value
#'
#'   For new observations selected from the additional dataset (`dataset_add`),
#'   `exist_flag` is set to the specified value.
#'
#' @param false_value False value
#'
#'   For new observations not selected from the additional dataset
#'   (`dataset_add`), `exist_flag` is set to the specified value.
#'
#'
#' @param filter Filter for observations to consider
#'
#'   *Deprecated*, please use the above `filter_add` argument instead.
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for selecting the first or last observation. If the argument is not
#'   specified, all observations are considered.
#'
#'   *Permitted Values*: a condition
#'
#' @inheritParams filter_extreme
#' @inheritParams derive_summary_records
#'
#' @details
#'   1. The additional dataset (`dataset_add`) is restricted as specified by the
#'   `filter_add` argument.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` argument) the first or last observation (with respect to the
#'   order specified for the `order` argument and the mode specified for the
#'   `mode` argument) is selected.
#'   1. If `dataset_ref` is specified, observations which are in `dataset_ref`
#'   but not in the selected records are added.
#'   1. The variables specified by the `set_values_to` argument are added to
#'   the selected observations.
#'   1. The observations are added to input dataset.
#'
#'
#' @return The input dataset with the first or last observation of each by group
#'   added as new observations.
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
#' adlb <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
#'   "1",      1,          113,      1,
#'   "1",      2,          113,      2,
#'   "1",      3,          117,      3,
#'   "2",      1,          101,      1,
#'   "2",      2,          101,      2,
#'   "2",      3,           95,      3
#' )
#'
#' # Add a new record for each USUBJID storing the minimum value (first AVAL).
#' # If multiple records meet the minimum criterion, take the first value by
#' # AVISITN. Set AVISITN = 97 and DTYPE = MINIMUM for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL, AVISITN),
#'   mode = "first",
#'   filter_add = !is.na(AVAL),
#'   set_values_to = exprs(
#'     AVISITN = 97,
#'     DTYPE = "MINIMUM"
#'   )
#' )
#'
#' # Add a new record for each USUBJID storing the maximum value (last AVAL).
#' # If multiple records meet the maximum criterion, take the first value by
#' # AVISITN. Set AVISITN = 98 and DTYPE = MAXIMUM for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(desc(AVAL), AVISITN),
#'   mode = "first",
#'   filter_add = !is.na(AVAL),
#'   set_values_to = exprs(
#'     AVISITN = 98,
#'     DTYPE = "MAXIMUM"
#'   )
#' )
#'
#' # Add a new record for each USUBJID storing for the last value.
#' # Set AVISITN = 99 and DTYPE = LOV for these new records.
#' derive_extreme_records(
#'   adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISITN = 99,
#'     DTYPE = "LOV"
#'   )
#' )
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
#' derive_extreme_records(
#'   adrs,
#'   dataset_ref = adsl,
#'   dataset_add = adrs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   filter_add = PARAMCD == "OVR" & AVALC == "PD",
#'   order = exprs(ADT),
#'   exist_flag = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   mode = "first",
#'   set_values_to = exprs(
#'     PARAMCD = "PD",
#'     PARAM = "Disease Progression",
#'     AVAL = yn_to_numeric(AVALC),
#'     ANL01FL = "Y",
#'     ADT = ADT
#'   )
#' )
#'
#' # derive parameter indicating death
#' derive_extreme_records(
#'   dataset_ref = adsl,
#'   dataset_add = adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   filter_add = !is.na(DTHDT),
#'   exist_flag = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   mode = "first",
#'   set_values_to = exprs(
#'     PARAMCD = "DEATH",
#'     PARAM = "Death",
#'     ANL01FL = "Y",
#'     ADT = DTHDT
#'   )
#' )
derive_extreme_records <- function(dataset = NULL,
                                   dataset_add = NULL,
                                   dataset_ref = NULL,
                                   by_vars = NULL,
                                   order = NULL,
                                   mode = NULL,
                                   filter_add = NULL,
                                   check_type = "warning",
                                   exist_flag = NULL,
                                   true_value = "Y",
                                   false_value = "N",
                                   set_values_to,
                                   filter) {
  if (!missing(filter)) {
    deprecate_warn(
      "0.11.0",
      "derive_extreme_records(filter = )",
      "derive_extreme_records(filter_add = )"
    )
    filter_add <- enexpr(filter)
  }

  # Check input arguments
  assert_vars(by_vars, optional = is.null(dataset_ref))
  assert_expr_list(order, optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = expr_c(
      by_vars, extract_vars(order)
    ),
    optional = TRUE
  )
  mode <- assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE,
    optional = TRUE
  )
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  exist_flag <- assert_symbol(enexpr(exist_flag), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_varval_list(set_values_to)
  if (is.null(dataset) && is.null(dataset_add)) {
    abort(paste(
      "Neither `dataset` nor `dataset_add` is specified.",
      "At least one of them must be specified.",
      sep = "\n"
    ))
  }

  # Create new observations
  if (is.null(dataset_add)) {
    dataset_add <- dataset
  }
  new_add_obs <- filter_if(dataset_add, filter_add)

  if (!is.null(order)) {
    new_add_obs <- filter_extreme(
      new_add_obs,
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    )
  }

  if (!is.null(dataset_ref)) {
    add_vars <- colnames(dataset_add)
    ref_vars <- colnames(dataset_ref)

    new_ref_obs <- anti_join(
      select(dataset_ref, intersect(add_vars, ref_vars)),
      select(new_add_obs, !!!by_vars),
      by = map_chr(by_vars, as_name)
    )

    if (!is.null(exist_flag)) {
      new_add_obs <- mutate(new_add_obs, !!exist_flag := true_value)
      new_ref_obs <- mutate(new_ref_obs, !!exist_flag := false_value)
    }
    new_obs <- bind_rows(new_add_obs, new_ref_obs)
  } else {
    new_obs <- new_add_obs
  }

  new_obs <- process_set_values_to(
    new_obs,
    set_values_to = set_values_to
  )

  # Create output dataset
  bind_rows(dataset, new_obs)
}
