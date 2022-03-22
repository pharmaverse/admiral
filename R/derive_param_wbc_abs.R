#' Add a parameter for lab differentials converted to absolute values
#'
#' @description
#'
#' Add a parameter by converting lab differentials from percentage to absolute values
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` argument are expected to be present.
#'
#'   The variable specified by `by_vars` and `LBTESTCD` must be a unique key of
#'   the input dataset after restricting it by the `filter_diff` condition,
#'   and to the parameters specified by `wbc_code` and `diff_code`.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `vars()` defining the variables to be set for the
#'   new parameter, e.g. `vars(PARAMCD = "LYMPHSI",
#'   PARAM = "Lymphocytes Abs (10^9/L) (SI Result)")` is
#'   expected.
#'
#' @param filter_diff Filter condition to identify differential records
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, only the differential records specified by `diff_code` and satisfy the
#'   `filter_diff` condition are kept.
#'
#'   Default: `LBSTRESU == "fraction of 1"`
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param wbc_code White Blood Cell (WBC) parameter
#'
#'   The observations where `LBTESTCD` equals the specified value are considered
#'   as the WBC absolute results to use for converting the differentials.
#'
#'   Default: `"WBC"`
#'
#'   Permitted Values: character value
#'
#' @param diff_code white blood differential parameter
#'
#'   The observations where `LBTESTCD` equals the specified value are considered
#'   as the white blood differential lab results in percentage value (i.e. `LBSTRESU`
#'   is "fraction of 1") to be converted into absolute value.
#'
#'
#' @details
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{White Blood Cell Count  * Percentage Value}{100}}.
#' New records are created for each group of records (grouped by `by_vars`) if 1) the white blood
#' cell component in absolute value is not already available from the input dataset, and 2) the
#' white blood cell absolute value (identified by `wbc_code`) and the white blood cell differential
#' (identified by `diff_code`) are both present.
#'
#' @author Annie Yang
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation bds adlb
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' test_lb <- tibble::tribble(
#' ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
#' "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
#' "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
#' "P01",    "LYMLE",   90,   "fraction of 1", "CYCLE 1 DAY 1",
#' "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1",
#' "P01",    "ALB",     36,   "g/L",           "CYCLE 2 DAY 1",
#' "P02",    "WBC",     33,   "10^9/L",        "CYCLE 1 DAY 1",
#' "P02",    "LYMPHSI", 29,   "10^9/L",        "CYCLE 1 DAY 1",
#' "P02",    "LYMLE",   87,   "fraction of 1", "CYCLE 1 DAY 1",
#' "P03",    "LYMLE",   89,   "fraction of 1", "CYCLE 1 DAY 1"
#' )
#'
#' derive_param_wbc_abs (dataset = test_lb,
#'                       by_vars = vars(USUBJID, VISIT),
#'                       set_values_to = vars(PARAMCD = "LYMPHSI",
#'                                            PARAM = "Lymphocytes Abs (10^9/L) (SI Result)",
#'                                            DTYPE = "CALCULATION",
#'                                            AVALU = "10^9/L"),
#'                       filter_diff = LBSTRESU == "fraction of 1",
#'                       wbc_code = "WBC",
#'                       diff_code = "LYMLE")

derive_param_wbc_abs <- function(dataset,
                                 by_vars,
                                 set_values_to,
                                 filter_diff = LBSTRESU == "fraction of 1",
                                 wbc_code = "WBC",
                                 diff_code) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, LBTESTCD, LBSTRESN, LBSTRESU))
  assert_character_scalar(wbc_code)
  assert_character_scalar(diff_code)
  filter_diff <- assert_filter_cond(enquo(filter_diff))
  by_vars_str <- vars2chr(by_vars)

  # select observations and variables required for new observations.
  dataset_temp <- dataset %>%
    filter((LBTESTCD == !!wbc_code) |
             (LBTESTCD == !!diff_code & !!filter_diff) |
             LBTESTCD == !!quo_get_expr(set_values_to$PARAMCD)) %>%
    select(by_vars_str, LBTESTCD, LBSTRESN, LBSTRESU)

  # only keep records where absolute value do not already exist.
  # add PARAMCD and AVAL (required in input dataset for derive_derived_param() function).
  dataset_abs <- dataset_temp %>%
    filter(LBTESTCD == !!quo_get_expr(set_values_to$PARAMCD)) %>%
    mutate(temp_flag = "Y") %>%
    select(by_vars_str, temp_flag)

  dataset_temp <- dataset_temp %>%
    left_join(dataset_abs, by = by_vars_str) %>%
    filter(is.na(temp_flag)) %>%
    mutate(PARAMCD = LBTESTCD, AVAL = LBSTRESN)

  # create new parameters.
  temp_set_values_to <- quo_c(set_values_to, vars(temp_var = "new"))

  dataset_new <- dataset_temp %>%
    derive_derived_param(
      parameters = c(wbc_code, diff_code),
      by_vars = by_vars,
      analysis_value = ((!!sym(paste0("AVAL.", wbc_code)) *
                           !!sym(paste0("AVAL.", diff_code))) / 100),
      set_values_to = temp_set_values_to
    ) %>%
    filter(temp_var == "new") %>%
    select(-starts_with("temp_"))

  # if no new records are added, output note and return original dataset,
  # else append new records to the original input dataset.
  if (nrow(dataset_new) == 0L) {
    message("No source records meet condition for calculation, therefore no new records created")
    dataset
  } else {
    bind_rows(dataset, dataset_new)
  }

}
