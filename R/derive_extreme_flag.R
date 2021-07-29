#' Adds a Variable Flagging the First or Last Observation Within Each By Group
#'
#' Adds a variable flagging the first or last observation within each by group
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_var Variable to add
#'
#'   The specified variable is added to the output dataset. It is set to `"Y"`
#'   for the first or last observation (depending on the mode) of each by group.
#'
#'   Permitted Values: list of name-value pairs
#'
#' @param order Sort order
#'
#'   The first or last observation is determined with respect to the specified
#'   order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param mode Flag mode
#'
#'   Determines of the first or last observation is flagged.
#'
#'   Permitted Values: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param flag_filter Filter for flag data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the parameter is not specified, all observations are
#'   considered.
#'
#'   Permitted Values: a condition
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   Default: `"warning"`
#'
#'   Permitted Values: `"none"`, `"warning"`, `"error"`
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the flag mode) is included in
#'   the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("vs")
#'
#' # Flag last value for each patient, test, and visit, baseline observations are ignored
#' vs %>%
#'   derive_extreme_flag(
#'     new_var = LASTFL,
#'     by_vars = vars(USUBJID, VSTESTCD, VISIT),
#'     order = vars(VSTPTNUM),
#'     mode = "last",
#'     flag_filter = VISIT != "BASELINE"
#'   ) %>%
#'   arrange(USUBJID, VSTESTCD, VISITNUM, VSTPTNUM) %>%
#'   select(USUBJID, VSTESTCD, VISIT, VSTPTNUM, VSSTRESN, LASTFL)
#'
#' # Baseline (ABLFL) examples:
#'
#' input <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~DTYPE,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
#'   "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
#'   "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,
#'
#'   "TEST01", "PAT02",  "PARAM01", "SCREENING",as.Date("2021-04-27"), 15.0, "AVERAGE",
#'   "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE",
#'   "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
#'   "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
#'   "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE",
#'
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, "AVERAGE",
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-25"), 14.0, "AVERAGE",
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-23"), 15.0, NA,
#'   "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "AVERAGE",
#'   "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,
#'
#'   "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, NA,
#'   "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
#'   "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
#'   "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
#'   "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA
#' )
#'
#' # Last observation
#' derive_extreme_flag(
#'   input,
#'   new_var = ABLFL,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   order = vars(ADT),
#'   mode = "last",
#'   flag_filter = AVISIT == "BASELINE"
#' )
#'
#' # Worst observation - Direction = High
#' derive_extreme_flag(
#'   input,
#'   new_var = ABLFL,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   order = vars(AVAL, ADT),
#'   mode = "last",
#'   flag_filter = AVISIT == "BASELINE"
#' )
#'
#' # Worst observation - Direction = Lo
#' derive_extreme_flag(
#'   input,
#'   new_var = ABLFL,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   order = vars(desc(AVAL), ADT),
#'   mode = "last",
#'   flag_filter = AVISIT == "BASELINE"
#' )
#'
#' # Average observation
#' derive_extreme_flag(
#'   input,
#'   new_var = ABLFL,
#'   by_vars = vars(USUBJID, PARAMCD),
#'   order = vars(ADT, desc(AVAL)),
#'   mode = "last",
#'   flag_filter = AVISIT == "BASELINE" & DTYPE == "AVERAGE"
#' )
#'
derive_extreme_flag <- function(dataset,
                                new_var,
                                by_vars,
                                order,
                                mode,
                                flag_filter = NULL,
                                check_type = "warning") {
  # checking and quoting
  new_var <- assert_symbol(enquo(new_var))
  assert_vars(by_vars)
  assert_order_vars(order)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, !!!extract_vars(order)))
  assert_character_scalar(mode, values = c("first", "last"))
  flag_filter <- assert_filter_cond(enquo(flag_filter), optional = TRUE)
  assert_character_scalar(check_type, values = c("none", "warning", "error"))

  # select data to consider for flagging
  if (!quo_is_null(flag_filter)) {
    data <- dataset %>% filter(!!flag_filter)
    data_ignore <- dataset %>%
      filter(!(!!flag_filter) | is.na(!!flag_filter))
  } else {
    data <- dataset
  }

  # create flag
  data <- data %>%
    derive_obs_number(new_var = temp_obs_nr,
                      order = order,
                      by_vars = by_vars,
                      check_type = check_type)

  if (mode == "first") {
    data <- data %>%
      mutate(!!new_var := if_else(temp_obs_nr == 1, "Y", NA_character_))
  } else {
    data <- data %>%
      group_by(!!!by_vars) %>%
      mutate(!!new_var := if_else(temp_obs_nr == n(), "Y", NA_character_)) %>%
      ungroup()
  }

  # add ignored data
  if (!quo_is_null(flag_filter)) {
    data <- data %>% bind_rows(data_ignore)
  }

  # remove temporary variable
  data %>% select(-temp_obs_nr)
}

#' Adds a Variable Flagging the maximal / minimal value within a group of observations
#'
#' @inheritParams derive_extreme_flag
#' @param dataset Input dataset.
#' Variables specified by `by_vars`, `order`, `param_var`, and `analysis_var` are expected.
#' @param new_var Variable to add to the `dataset`.
#' It is set `"Y"` for the maximal / minimal observation of each group,
#' see Details section for more information.
#' @param order Sort order.
#' Used to determine maximal / minimal observation if they are not unique,
#' see Details section for more information.
#' @param param_var Variable with the parameter values for which the maximal / minimal
#' value is calculated.
#' @param analysis_var Variable with the measurement values for which the maximal / minimal
#' value is calculated.
#' @param worst_high Character with `param_var` values specifying the parameters
#' referring to "high".
#' Use `character(0)` if not required.
#' @param worst_low Character with `param_var` values specifying the parameters
#' referring to "low".
#' Use `character(0)` if not required.
#'
#' @details For each group with respect to the variables specified by the `by_vars` parameter,
#' the maximal / minimal observation of `analysis_var`
#' is labelled in the `new_var` column as `"Y"`
#' if its `param_var` is in `worst_high` / `worst_low`,
#' otherwise it is assigned `NA`.
#' If there is more than one such maximal / minimal observation,
#' the first one with respect to the order specified by the `order` parameter is flagged.
#'
#' @author Ondrej Slama
#'
#' @return The input dataset with the new flag variable added.
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#'
#' input <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
#'   "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
#'   "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
#'   "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,
#'
#'   "TEST01", "PAT02",  "PARAM01", "SCREENING",as.Date("2021-04-27"), 15.0,
#'   "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
#'   "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
#'   "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
#'   "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,
#'
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-25"), 14.0,
#'   "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-23"), 15.0,
#'   "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
#'   "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,
#'
#'   "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
#'   "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
#'   "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
#'   "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
#'   "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,
#'
#'   "TEST01", "PAT02",  "PARAM03", "SCREENING",as.Date("2021-04-27"), 15.0,
#'   "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-25"), 14.0,
#'   "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-23"), 15.0,
#'   "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-27"), 10.0,
#'   "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
#' )
#'
#' derive_worst_flag(
#'   input,
#'   new_var = WORSTFL,
#'   by_vars = vars(USUBJID, PARAMCD, AVISIT),
#'   order = vars(desc(ADT)),
#'   param_var = PARAMCD,
#'   analysis_var = AVAL,
#'   worst_high = c("PARAM01", "PARAM03"),
#'   worst_low = "PARAM02"
#' )
#'
#'\dontrun{
#' # example with ADVS
#' derive_worst_flag(
#'   advs,
#'   new_var = WORSTFL,
#'   by_vars = vars(USUBJID, PARAMCD, AVISIT),
#'   order = vars(ADT, ATPTN),
#'   param_var = PARAMCD,
#'   analysis_var = AVAL,
#'   worst_high = c("SYSBP", "DIABP"),
#'   worst_low = "RESP",
#'   flag_filter = !is.na(AVISIT) & !is.na(AVAL)
#' )
#'}
#'
derive_worst_flag <- function(dataset,
                              new_var,
                              by_vars,
                              order,
                              param_var,
                              analysis_var,
                              worst_high,
                              worst_low,
                              flag_filter = NULL,
                              check_type = "warning") {

  # perform argument checks
  new_var <- assert_symbol(enquo(new_var))
  param_var <- assert_symbol(enquo(param_var))
  analysis_var <- assert_symbol(enquo(analysis_var))
  assert_vars(by_vars)
  assert_order_vars(order)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, extract_vars(order), param_var, analysis_var)
  )
  assert_character_vector(worst_high)
  assert_character_vector(worst_low)
  flag_filter <- assert_filter_cond(enquo(flag_filter), optional = TRUE)

  # additional checks for worstflag - parameters overlap
  if (length(intersect(worst_high, worst_low)) > 0) {
    err_msg <- paste(
      "The following parameter(-s) are both assigned to `worst_high` and `worst_low` flags:",
      paste0(intersect(worst_high, worst_low), collapse = ", ")
    )
    abort(err_msg)
  }

  # additional checks for worstflag - parameters not available
  param_var_str <- as_string(quo_get_expr(param_var))
  if (length(worst_high) > 0 &&
      !all(worst_high %in% dataset[[param_var_str]])) {
    err_msg <- paste0(
      "The following parameter(-s) in `worst_high` are not available in column ",
      param_var_str,
      ": ",
      paste0(worst_high[!worst_high %in% dataset[[param_var_str]]], collapse = ", ")
    )
    abort(err_msg)
  }

  # additional checks for worstflag - parameters not available
  if (length(worst_low) > 0 &&
      !all(worst_low %in% dataset[[param_var_str]])) {
    err_msg <- paste0(
      "The following parameter(-s) in `worst_low` are not available in column ",
      param_var_str,
      ": ",
      paste0(worst_low[!worst_low %in% dataset[[param_var_str]]], collapse = ", ")
    )
    abort(err_msg)
  }

  # derive worst-flag
  bind_rows(
    derive_extreme_flag(
      dataset = filter(dataset, !!param_var %in% worst_low),
      new_var = !!new_var,
      by_vars = by_vars,
      order = quo_c(analysis_var, order),
      mode = "first",
      flag_filter = !!flag_filter,
      check_type = check_type
    ),
    derive_extreme_flag(
      dataset = filter(dataset, !!param_var %in% worst_high),
      new_var = !!new_var,
      by_vars = by_vars,
      order = quo_c(quo(desc(!!quo_get_expr(analysis_var))), order),
      mode = "first",
      flag_filter = !!flag_filter,
      check_type = check_type
    ),
    filter(dataset, !(!!param_var %in% c(worst_low, worst_high)))
  )
}
