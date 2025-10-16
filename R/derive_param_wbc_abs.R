#' Add a parameter for lab differentials converted to absolute values
#'
#' @description
#'
#' Add a parameter by converting lab differentials from fraction or percentage to absolute values
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset, and to the parameters specified by `wbc_code` and `diff_code`.
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `exprs()` defining the variables to be set for the
#'   new parameter, e.g. `exprs(PARAMCD = "LYMPH", PARAM = "Lymphocytes Abs (10^9/L)")` is
#'   expected.
#'
#' @param get_unit_expr An expression providing the unit of the parameter
#'
#'   The result is used to check the units of the input parameters.
#'
#' @permitted a variable containing unit from the input dataset, or a function call,
#'   for example, `get_unit_expr = extract_unit(PARAM)`.
#'
#' @param wbc_unit A string containing the required unit of the WBC parameter
#'
#'
#' @param diff_type A string specifying the type of differential
#'
#' @permitted `"percent"`, `"fraction"`
#'
#' @param by_vars Grouping variables
#'
#' `r roxygen_param_by_vars()`
#'
#' @param wbc_code White Blood Cell (WBC) parameter
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the WBC absolute results to use for converting the differentials.
#'
#'
#' @permitted character value
#'
#' @param diff_code white blood differential parameter
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the white blood differential lab results in fraction or percentage value to be converted
#'   into absolute value.
#'
#'
#' @details
#' If `diff_type` is `"percent"`, the analysis value of the new parameter is derived as
#' \deqn{\frac{White Blood Cell Count  * Percentage Value}{100}}
#'
#' If `diff_type` is `"fraction"`, the analysis value of the new parameter is derived as
#' \deqn{White Blood Cell Count  * Fraction Value}
#'
#'
#' New records are created for each group of records (grouped by `by_vars`) if 1) the white blood
#' cell component in absolute value is not already available from the input dataset, and 2) the
#' white blood cell absolute value (identified by `wbc_code`) and the white blood cell differential
#' (identified by `diff_code`) are both present.
#'
#'
#' @return The input dataset with the new parameter added
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' test_lb <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
#'   "P01", "WBC", 33, "Leukocyte Count (10^9/L)", "CYCLE 1 DAY 1",
#'   "P01", "WBC", 38, "Leukocyte Count (10^9/L)", "CYCLE 2 DAY 1",
#'   "P01", "LYMLE", 0.90, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
#'   "P01", "LYMLE", 0.70, "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
#'   "P01", "ALB", 36, "Albumin (g/dL)", "CYCLE 2 DAY 1",
#'   "P02", "WBC", 33, "Leukocyte Count (10^9/L)", "CYCLE 1 DAY 1",
#'   "P02", "LYMPH", 29, "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1",
#'   "P02", "LYMLE", 0.87, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
#'   "P03", "LYMLE", 0.89, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1"
#' )
#'
#' derive_param_wbc_abs(
#'   dataset = test_lb,
#'   by_vars = exprs(USUBJID, VISIT),
#'   set_values_to = exprs(
#'     PARAMCD = "LYMPH",
#'     PARAM = "Lymphocytes Abs (10^9/L)",
#'     DTYPE = "CALCULATION"
#'   ),
#'   get_unit_expr = extract_unit(PARAM),
#'   wbc_code = "WBC",
#'   diff_code = "LYMLE",
#'   diff_type = "fraction"
#' )
derive_param_wbc_abs <- function(dataset,
                                 by_vars,
                                 set_values_to,
                                 get_unit_expr,
                                 wbc_unit = "10^9/L",
                                 wbc_code = "WBC",
                                 diff_code,
                                 diff_type = "fraction") {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, PARAMCD, AVAL))
  assert_character_scalar(wbc_code)
  assert_character_scalar(diff_code)
  assert_character_scalar(diff_type, values = c("percent", "fraction"))
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  assert_character_scalar(wbc_unit)
  by_vars_str <- vars2chr(by_vars)

  # Check unit of input parameters
  assert_unit(dataset, wbc_code, required_unit = wbc_unit, get_unit_expr = !!get_unit_expr)

  # Select observations and variables required for new observations.
  dataset_temp <- dataset %>%
    filter(
      PARAMCD == !!wbc_code |
        PARAMCD == !!diff_code |
        PARAMCD == !!set_values_to$PARAMCD
    ) %>%
    select(!!!by_vars, PARAMCD, AVAL)

  # Only keep records where absolute value do not already exist.
  temp_flag <- get_new_tmp_var(dataset_temp, prefix = "tmp_flag")
  dataset_abs <- dataset_temp %>%
    filter(PARAMCD == !!set_values_to$PARAMCD) %>%
    mutate(!!temp_flag := "Y") %>%
    select(!!!by_vars, !!temp_flag)

  dataset_temp <- dataset_temp %>%
    left_join(dataset_abs, by = by_vars_str) %>%
    filter(is.na(!!temp_flag)) %>%
    remove_tmp_vars()

  if (diff_type == "percent") {
    analysis_value <- expr(!!sym(paste0("AVAL.", wbc_code)) *
      !!sym(paste0("AVAL.", diff_code)) / 100)
  } else if (diff_type == "fraction") {
    analysis_value <- expr(!!sym(paste0("AVAL.", wbc_code)) *
      !!sym(paste0("AVAL.", diff_code)))
  }

  # Create new parameter.
  withCallingHandlers(
    dataset_new <- dataset_temp %>%
      derive_param_computed(
        parameters = c(
          wbc_code,
          diff_code
        ),
        by_vars = by_vars,
        set_values_to = exprs(
          AVAL = !!analysis_value,
          !!!set_values_to
        )
      ) %>%
      filter(PARAMCD == !!set_values_to$PARAMCD),
    derive_param_computed_all_na = function(cnd) {
      cnd_muffle(cnd)
    }
  )

  # If no new records are added, output note and return original dataset,
  # else append new records to the original input dataset.
  if (nrow(dataset_new) == 0L) {
    cli_inform("No source records meet condition for calculation, therefore no new records created")
    dataset
  } else {
    bind_rows(dataset, dataset_new)
  }
}
