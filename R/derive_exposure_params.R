#' Add a Record Computed from the Aggregated Analysis Value of another Parameter
#'
#' Add a record computed from the aggregated analysis value of another parameter.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, `AVAL`, `AVALC`, `ASTDTM`
#'   and `AENDTM` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `parameters`.
#'
#' @param filter Filter condition
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param input_param Required parameter code
#'
#' The observations where `PARAMCD` equals the specified value are considered to compute the
#' summary record.
#'
#'   *Permitted Values:* A character of `PARAMCD` value
#'
#' @param fns List of formulas specifying variable to use for aggregations.
#'
#' This can include base functions like `mean()`, `min()`, `max()`, `median()`,
#' `sd()`, or `sum()` or any other user-defined aggregation function.
#' For example, `fns = list(AVAL~ mean)`.
#'
#'   In general,
#'
#'   + LHS refer to the one or more variable to use for summarizing.
#'   + RHS refer to a **single** summary function.
#'
#'   In the formula representation e.g., `CHG ~ sum(., na.rm = TRUE)`, a `.`
#'   serves as the data to be summarized which refers to the variable `CHG`.
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset.
#'
#'   *Permitted Values:* list of variables
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations (e.g. `vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL")`).
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @author Samia Kabi
#'
#' @return The input dataset with a new record added. ASTDT(M) and AENDT(M) are computed as the
#' minimum/maximum dates by `by_vars`
#'
#' @keywords derivation bds
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(lubridate)
#' library(stringr)
#' adex <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC, ~VISIT, ~ASTDT, ~AENDT,
#'   "01-701-1015", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "01-701-1015", "DOSE", 85, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "01-701-1015", "DOSE", 82, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "01-701-1015", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "01-701-1015", "ADJ", NA, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "01-701-1015", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "01-701-1017", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "01-701-1017", "DOSE", 50, NA_character_, "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "01-701-1017", "DOSE", 65, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02"),
#'   "01-701-1017", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "01-701-1017", "ADJ", NA, "ADVERSE EVENT", "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "01-701-1017", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02")
#' ) %>%
#'   mutate(ASTDTM = ymd_hms(paste(ASTDT, "00:00:00")), AENDTM = ymd_hms(paste(AENDT, "00:00:00")))
#'
#'   # Cumulative dose
#'   adex %>%
#'   derive_exposure_params(
#'     by_vars = vars(USUBJID),
#'     set_values_to = vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL"),
#'     input_param = "DOSE",
#'     fns = AVAL ~ sum(., na.rm = TRUE)
#'   ) %>%
#'   # average dose in w2-24
#'   derive_exposure_params(
#'     by_vars = vars(USUBJID),
#'     filter = VISIT %in% c("WEEK2", "WEEK 24"),
#'     set_values_to = vars(PARAMCD = "AVDW224", PARCAT1 = "WEEK2-24"),
#'     input_param = "DOSE",
#'     fns = AVAL ~ mean(., na.rm = TRUE)
#'   ) %>%
#'   # Any dose adjustement?
#'   derive_exposure_params(
#'     by_vars = vars(USUBJID),
#'     set_values_to = vars(PARAMCD = "TADJ", PARCAT1 = "OVERALL"),
#'     input_param = "ADJ",
#'     fns = AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)
#'   )



derive_exposure_params <- function(dataset,
                                   by_vars,
                                   input_param,
                                   fns,
                                   filter = NULL,
                                   set_values_to = NULL) {
  by_vars <- assert_vars(by_vars)
  assert_data_frame(dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL, AVALC, ASTDTM, AENDTM))
  )
  assert_character_scalar(input_param)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(input_param, values = params_available)

  subset_ds <- dataset %>%
    filter_if(filter)

  add_data <- subset_ds %>%
    filter(PARAMCD == input_param) %>%
    derive_summary_records(
      by_vars = by_vars,
      fns = list(fns),
      set_values_to = set_values_to
    ) %>%
    filter(PARAMCD == quo_get_expr(set_values_to$PARAMCD))

  # add the dates for the derived parameters
  by_vars <- vars2chr(by_vars)
  dates <- subset_ds %>%
    group_by(!!!syms(by_vars)) %>%
    summarise(
      ASTDTM___ = min(ASTDTM, na.rm = TRUE),
      AENDTM___ = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
    )

  expo_data <- add_data %>%
    left_join(dates, by = by_vars) %>%
    mutate(
      ASTDTM = coalesce(as_iso_dttm(ASTDTM), as_iso_dttm(ASTDTM___)),
      AENDTM = coalesce(as_iso_dttm(AENDTM), as_iso_dttm(AENDTM___)),
      ASTDT = date(ASTDTM),
      AENDT = date(AENDTM),
      !!!set_values_to
    ) %>%
    select(-ends_with("___"))

  data <- bind_rows(dataset, expo_data)
}
