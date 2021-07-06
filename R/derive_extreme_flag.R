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
