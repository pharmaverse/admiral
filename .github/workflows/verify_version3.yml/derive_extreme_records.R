#' Add the First or Last Observation for Each By Group as New Records
#'
#' Add the first or last observation for each by group as new observations. The
#' new observations can be selected from the additional dataset. This function can
#' be used for adding the maximum or minimum value as a separate visit.
#' All variables of the selected observation are kept. This distinguishes
#' `derive_extreme_records()` from `derive_summary_records()`,
#' where only the by variables are populated for the new records.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' If the argument is not specified (or set to `NULL`), a new dataset is
#' created. Otherwise, the new records are appended to the specified dataset.
#'
#' @permitted [dataset]
#'
#' @param dataset_ref Reference dataset
#'
#'   The variables specified for `by_vars` are expected. For each observation of
#'   the specified dataset a new observation is added to the input dataset.
#'
#'   For records which are added from `dataset_ref` because there are no records
#'   in `dataset_add` for the by group only those variables are kept which are
#'   also in `dataset_add` (and are included in `keep_source_vars`).
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The additional dataset, which determines the by groups returned in the input dataset,
#'   based on the groups that exist in this dataset after being subset by `filter_add`.
#'
#'   The variables specified in the `by_vars` and `filter_add` parameters are expected
#'   in this dataset. If `mode` and `order` are specified, the first or last observation
#'   within each by group, defined by `by_vars`, is selected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   If `dataset_ref` is specified, this argument must be specified.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations in `dataset_add` fulfilling the specified condition are
#'   considered.
#'
#' @permitted [condition]
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is added
#'   to the input dataset. If `"last"` is specified, the last observation of
#'   each by group is added to the input dataset.
#'
#' @permitted [mode]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) additional dataset are not unique
#'   with respect to the by variables and the order.
#'
#' @permitted [msg_type]
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
#' @permitted Variable name
#'
#' @param true_value True value
#'
#'   For new observations selected from the additional dataset (`dataset_add`),
#'   `exist_flag` is set to the specified value.
#'
#' @permitted [char_scalar]
#'
#' @param false_value False value
#'
#'   For new observations not selected from the additional dataset
#'   (`dataset_add`), `exist_flag` is set to the specified value.
#'
#' @permitted [char_scalar]
#'
#' @param keep_source_vars Variables to be kept in the new records
#'
#'   A named list or tidyselect expressions created by `exprs()` defining the
#'   variables to be kept for the new records. The variables specified for
#'   `by_vars` and `set_values_to` need not be specified here as they are kept
#'   automatically.
#'
#' @permitted [var_list_tidyselect]
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
#'   but not in the selected records are added. Variables that are common
#'   across `dataset_ref`, `dataset_add` and `keep_source_vars()` are also
#'   populated for the new observations.
#'   1. The variables specified by the `set_values_to` argument are added to
#'   the selected observations.
#'   1. The variables specified by the `keep_source_vars` argument are selected
#'   along with the variables specified in `by_vars` and `set_values_to`
#'   arguments.
#'   1. The observations are added to input dataset (`dataset`). If no input
#'   dataset is provided, a new dataset is created.
#'
#'
#' @return The input dataset with the first or last observation of each by group
#'   added as new observations.
#'
#' @seealso [derive_summary_records()]
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examplesx
#'
#' @caption Add last/first record as new record
#' @info For each subject the last record should be added as a new visit.
#'
#' - The source dataset for the new records is specified by the `dataset_add`
#'   argument. Here it is the same as the input dataset.
#' - The groups for which new records are added are specified by the `by_vars`
#'   argument. Here for each *subject* a record should be added. Thus `by_vars =
#'   exprs(USUBJID)` is specified.
#' - As there are multiple records per subject, the `mode` and `order` arguments
#'   are specified to request that the *last* record is selected when the
#'   records are sorted by visit (`AVISITN`). The records are sorted by each by
#'   group (`by_vars`) separately, i.e., it is not necessary to add the
#'   variables from `by_vars` to `order`.
#' - To avoid duplicates in the output dataset the `set_values_to` argument is
#'   specified to set the visit (`AVISIT`) to a special value for the *new*
#'   records.
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' adlb <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVAL,
#'   "1",      1,          113,
#'   "1",      2,          111,
#'   "2",      1,          101,
#'   "2",      2,           NA,
#'   "3",      1,           NA,
#' )
#'
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISITN = 99
#'   )
#' )
#'
#' @caption Restricting source records (`filter_add`)
#' @info The source records can be restricted by the `filter_add` argument,
#'   e.g., to exclude visits with missing analysis value from selecting for the
#'   last visit record:
#' @code
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISITN = 99
#'   )
#' )
#' @info Please note that new records are added only for subjects in the
#'   *restricted* source data. Therefore no new record is added for subject
#'   `"3"`.
#'
#' @caption Adding records for groups not in source (`dataset_ref`)
#' @info Adding records for groups which are not in the source data can be
#'   achieved by specifying a reference dataset by the `dataset_ref` argument.
#'   For example, specifying the input dataset for `dataset_ref` below ensures
#'   that new records are added also for subject without a valid analysis value:
#' @code
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL),
#'   dataset_ref = adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISITN = 99
#'   )
#' )
#'
#' @caption Selecting variables for new records (`keep_source_vars`)
#' @info Which variables from the source dataset are kept for the new records
#'   can be specified by the `keep_source_vars` argument. Variables specified by
#'   `by_vars` or `set_values_to` don't need to be added to `keep_source_vars`
#'   as these are always kept.
#' @code
#' adlb <- tribble(
#'   ~USUBJID, ~AVISIT,  ~AVAL, ~LBSEQ,
#'   "1",      "WEEK 1",   123,      1,
#'   "1",      "WEEK 2",   101,      2,
#'   "2",      "WEEK 1",    99,      1,
#'   "2",      "WEEK 2",   110,      2,
#'   "2",      "WEEK 3",    93,      3
#' )
#'
#' derive_extreme_records(
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL),
#'   mode = "first",
#'   keep_source_vars = exprs(AVAL),
#'   set_values_to = exprs(
#'     AVISIT = "MINIMUM"
#'   )
#' )
#'
#' @caption Handling duplicates (`check_type`)
#' @info The source records are checked regarding duplicates with respect to
#'   `by_vars` and `order`. By default, a warning is issued if any duplicates
#'   are found.
#' @code [expected_cnds = "duplicate_records"]
#' adlb <- tribble(
#'   ~USUBJID, ~AVISIT,  ~AVAL,
#'   "1",      "WEEK 1",   123,
#'   "1",      "WEEK 2",   123,
#'   "2",      "WEEK 1",    99,
#'   "2",      "WEEK 2",   110,
#'   "2",      "WEEK 3",    93,
#' )
#'
#' derive_extreme_records(
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL),
#'   mode = "first",
#'   set_values_to = exprs(
#'     AVISIT = "MINIMUM"
#'   )
#' )
#'
#' @info For investigating the issue, the dataset of the duplicate source records
#'   can be obtained by calling `get_duplicates_dataset()`:
#' @code
#' get_duplicates_dataset()
#' @info Common options to solve the issue are:
#' - Restricting the source records by specifying/updating the `filter_add` argument.
#' - Specifying additional variables for `order`.
#' - Setting `check_type = "none"` to ignore any duplicates.
#'
#' In this example it doesn't matter which of the records with the minimum value
#' is chosen because it doesn't affect the output dataset. Thus the third option
#' is used:
#' @code
#' derive_extreme_records(
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL),
#'   mode = "first",
#'   check_type = "none",
#'   set_values_to = exprs(
#'     AVISIT = "MINIMUM"
#'   )
#' )
#'
#' @caption Flagging existence of source records (`exist_flag`, `true_value`,
#'   `false_value`)
#' @info If the existence of a source record should be flagged, the `exist_flag`
#'   argument can be specified. The specified variable is set to `true_value` if
#'   a source record exists. Otherwise, it is set to `false_value`.
#'
#'   The `dataset_ref` argument should be specified as otherwise *all* new
#'   records originate from `dataset_add`, i.e., `exist_flag` would be set to
#'   `true_value` for all records.
#' @code
#' adsl <- tribble(
#'   ~USUBJID, ~DTHDT,
#'   "1",      ymd("2022-05-13"),
#'   "2",      ymd(""),
#'   "3",      ymd("")
#' )
#'
#' derive_extreme_records(
#'   dataset_ref = adsl,
#'   dataset_add = adsl,
#'   by_vars = exprs(USUBJID),
#'   filter_add = !is.na(DTHDT),
#'   exist_flag = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   set_values_to = exprs(
#'     PARAMCD = "DEATH",
#'     ADT = DTHDT
#'   )
#' )
#'
#' @caption Derive `DTYPE = "LOV"`
#' @info For each subject and parameter the last valid assessment (with respect
#'   to `AVISITN` and `LBSEQ`) should be selected and added as a new record to
#'   the input dataset. For the new records set `AVISIT = "PBL LAST"`, `AVISITN
#'   = 99`, and `DTYPE = "LOV"`.
#' @code
#' adlb <- tribble(
#'   ~USUBJID, ~AVISIT,    ~AVISITN, ~PARAMCD, ~AVAL, ~LBSEQ,
#'   "1",      "BASELINE",        1, "ABC",      120,      1,
#'   "1",      "WEEK 1",          2, "ABC",      113,      2,
#'   "1",      "WEEK 1",          2, "ABC",      117,      3,
#'   "2",      "BASELINE",        1, "ABC",      101,      1,
#'   "2",      "WEEK 1",          2, "ABC",      101,      2,
#'   "2",      "WEEK 2",          3, "ABC",       95,      3,
#'   "1",      "BASELINE",        1, "DEF",       17,      1,
#'   "1",      "WEEK 1",          2, "DEF",       NA,      2,
#'   "1",      "WEEK 1",          2, "DEF",       13,      3,
#'   "2",      "BASELINE",        1, "DEF",        9,      1,
#'   "2",      "WEEK 1",          2, "DEF",       10,      2,
#'   "2",      "WEEK 2",          3, "DEF",       12,      3
#' ) %>%
#' mutate(STUDYID = "XYZ", .before = USUBJID)
#'
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
#'   by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
#'   order = exprs(AVISITN, LBSEQ),
#'   mode = "last",
#'   set_values_to = exprs(
#'     AVISIT = "PBL LAST",
#'     AVISITN = 99,
#'     DTYPE = "LOV"
#'   )
#' )
#'
#' @caption Derive `DTYPE = "MINIMUM"`
#' @info For each subject and parameter the record with the minimum analysis
#'   value should be selected and added as a new record to the input dataset. If
#'   there are multiple records meeting the minimum, the first record with
#'   respect to `AVISIT` and `LBSEQ` should be selected. For the new records set
#'   `AVISIT = "PBL MIN"`, `AVISITN = 97`, and `DTYPE = "MINIMUM"`.
#' @code
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
#'   by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
#'   order = exprs(AVAL, AVISITN, LBSEQ),
#'   mode = "first",
#'   set_values_to = exprs(
#'     AVISIT = "PBL MIN",
#'     AVISITN = 97,
#'     DTYPE = "MINIMUM"
#'   )
#' )
#'
#' @caption Derive `DTYPE = "MAXIMUM"`
#' @info For each subject and parameter the record with the maximum analysis
#'   value should be selected and added as a new record to the input dataset. If
#'   there are multiple records meeting the maximum, the first record with
#'   respect to `AVISIT` and `LBSEQ` should be selected. For the new records set
#'   `AVISIT = "PBL MAX"`, `AVISITN = 98`, and `DTYPE = "MAXIMUM"`.
#' @code
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
#'   by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
#'   order = exprs(desc(AVAL), AVISITN, LBSEQ),
#'   mode = "first",
#'   set_values_to = exprs(
#'     AVISIT = "PBL MAX",
#'     AVISITN = 99,
#'     DTYPE = "MAXIMUM"
#'   )
#' )
#'
#' @caption Derive `DTYPE = "WOC"` or `DTYPE = "BOC"`
#' @info For each subject and parameter the record with the worst analysis
#'   value should be selected and added as a new record to the input dataset.
#'   The worst value is either the minimum or maximum value depending on the
#'   parameter. If there are multiple records meeting the worst value, the first
#'   record with respect to `AVISIT` and `LBSEQ` should be selected. For the new
#'   records set `AVISIT = "PBL WORST"`, `AVISITN = 96`, and `DTYPE = "WOC"`.
#'
#'   Here the maximum is considered worst for `PARAMCD = "ABC"` and the minimum
#'   for `PARAMCD = "DEF"`.
#' @code
#' derive_extreme_records(
#'   adlb,
#'   dataset_add = adlb,
#'   filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
#'   by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
#'   order = exprs(
#'     if_else(PARAMCD == "ABC", desc(AVAL), AVAL),
#'     AVISITN, LBSEQ
#'   ),
#'   mode = "first",
#'   set_values_to = exprs(
#'     AVISIT = "PBL WORST",
#'     AVISITN = 96,
#'     DTYPE = "WOC"
#'   )
#' )
#'
#' @caption Derive a parameter for the first disease progression (PD)
#' @info For each subject in the `ADSL` dataset a new parameter should be added
#'   to the input dataset which indicates whether disease progression (PD)
#'   occurred (set `AVALC = "Y"`, `AVAL = 1`) or not (set `AVALC = "N"`, `AVAL =
#'   0`). For the new parameter set `PARAMCD = "PD"` and `PARAM = "Disease
#'   Progression"`.
#' @code
#' adsl <- tribble(
#'   ~USUBJID, ~DTHDT,
#'   "1",      ymd("2022-05-13"),
#'   "2",      ymd(""),
#'   "3",      ymd("")
#' ) %>%
#'   mutate(STUDYID = "XX1234")
#'
#' adrs <- tribble(
#'   ~USUBJID, ~RSDTC,       ~AVALC, ~AVAL,
#'   "1",      "2020-01-02", "PR",       2,
#'   "1",      "2020-02-01", "CR",       1,
#'   "1",      "2020-03-01", "CR",       1,
#'   "2",      "2021-06-15", "SD",       3,
#'   "2",      "2021-07-16", "PD",       4,
#'   "2",      "2021-09-14", "PD",       4
#' ) %>%
#'   mutate(
#'     STUDYID = "XX1234", .before = USUBJID
#'   ) %>%
#'   mutate(
#'     ADT = ymd(RSDTC),
#'     PARAMCD = "OVR",
#'     PARAM = "Overall Response",
#'     .after = RSDTC
#'   )
#'
#' derive_extreme_records(
#'   adrs,
#'   dataset_ref = adsl,
#'   dataset_add = adrs,
#'   by_vars = get_admiral_option("subject_keys"),
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
#'   )
#' )
#'
#' @caption Derive parameter indicating death
#' @info For each subject in the `ADSL` dataset a new parameter should be
#'   created which indicates whether the subject died (set `AVALC = "Y"`, `AVAL
#'   = 1`) or not (set `AVALC = "N"`, `AVAL = 0`). For the new parameter set
#'   `PARAMCD = "DEATH"`, `PARAM = "Death"`, and `ADT` to the date of death
#'   (`DTHDT`).
#' @code
#' derive_extreme_records(
#'   dataset_ref = adsl,
#'   dataset_add = adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   filter_add = !is.na(DTHDT),
#'   exist_flag = AVALC,
#'   true_value = "Y",
#'   false_value = "N",
#'   mode = "first",
#'   keep_source_vars = exprs(AVALC),
#'   set_values_to = exprs(
#'     PARAMCD = "DEATH",
#'     PARAM = "Death",
#'     ADT = DTHDT
#'   )
#' )
#' @info The `keep_source_vars` argument is specified to avoid that all `ADSL`
#'   variables (like `DTHDT`) are copied to the parameter.
derive_extreme_records <- function(dataset = NULL,
                                   dataset_add,
                                   dataset_ref = NULL,
                                   by_vars = NULL,
                                   order = NULL,
                                   mode = NULL,
                                   filter_add = NULL,
                                   check_type = "warning",
                                   exist_flag = NULL,
                                   true_value = "Y",
                                   false_value = NA_character_,
                                   keep_source_vars = exprs(everything()),
                                   set_values_to) {
  # Check input arguments
  assert_vars(by_vars, optional = is.null(dataset_ref))
  assert_expr_list(order, optional = TRUE)
  assert_expr_list(keep_source_vars, optional = TRUE)

  assert_data_frame(
    dataset,
    optional = TRUE
  )
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(by_vars, extract_vars(order))
  )
  assert_data_frame(
    dataset_ref,
    required_vars = by_vars,
    optional = TRUE
  )
  mode <- assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE,
    optional = TRUE
  )
  exist_flag <- assert_symbol(enexpr(exist_flag), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_varval_list(set_values_to)

  # Create new observations
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

  new_obs <- new_obs %>%
    process_set_values_to(
      set_values_to = set_values_to
    ) %>%
    select(!!!by_vars, names(set_values_to), !!!keep_source_vars)


  # Create output dataset
  bind_rows(dataset, new_obs)
}
