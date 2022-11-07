#' Add an Aggregated Parameter and Derive the Associated Start and End Dates
#'
#' Add a record computed from the aggregated analysis value of another parameter and compute the
#' start (`ASTDT(M)`)and end date (`AENDT(M)`) as the minimum and maximum date by `by_vars`.
#'
#' @param dataset Input dataset
#'
#'   + The variables specified by the `by_vars`,`analysis_var` parameters and `PARAMCD` are
#'   expected,
#'   + Either `ASTDTM` and `AENDTM` or `ASTDT` and `AENDT` are also expected.
#'
#' @param filter Filter condition
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param input_code Required parameter code
#'
#' The observations where `PARAMCD` equals the specified value are considered to compute the
#' summary record.
#'
#'   *Permitted Values:* A character of `PARAMCD` value
#'
#' @param analysis_var Analysis variable.
#'
#' @param summary_fun Function that takes as an input the `analysis_var` and
#'   performs the calculation.
#'   This can include built-in functions as well as user defined functions,
#'   for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   *Permitted Values:* list of variables
#'
#' @param set_values_to Variable-value pairs
#'
#'   Set a list of variables to some specified value for the new observation(s)
#'   + LHS refer to a variable. It is expected that at least `PARAMCD` is defined.
#'   + RHS refers to the values to set to the variable. This can be a string, a symbol, a numeric
#'   value or NA.
#'   (e.g.  `vars(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")`).
#'   More general expression are not allowed.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @details For each group (with respect to the variables specified for the `by_vars` parameter),
#' an observation is added to the output dataset and the defined values are set to the defined
#' variables
#'
#'
#' @author Samia Kabi
#'
#' @return The input dataset with a new record added for each group (with respect to the variables
#' specified for the `by_vars` parameter). That is, a variable will only
#' be populated in this new record if it is specified in `by_vars`.
#' For each new record,
#' + the variable specified `analysis_var` is computed as defined by `summary_fun`,
#' + the variable(s) specified on the LHS of `set_values_to` are set to their paired value (RHS).
#' In addition, the start and end date are computed as the minimum/maximum dates by `by_vars`.
#'
#' If the input datasets contains
#' + both `AxxDTM` and `AxxDT` then all `ASTDTM`,`AENDTM`, `ASTDT`, `AENDT` are computed
#' + only `AxxDTM` then `ASTDTM`,`AENDTM` are computed
#' + only `AxxDT` then `ASTDT`,`AENDT` are computed.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#' library(stringr, warn.conflicts = FALSE)
#' adex <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC, ~VISIT, ~ASTDT, ~AENDT,
#'   "1015", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "1015", "DOSE", 85, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "1015", "DOSE", 82, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "1015", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
#'   "1015", "ADJ", NA, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
#'   "1015", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
#'   "1017", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "1017", "DOSE", 50, NA_character_, "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "1017", "DOSE", 65, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02"),
#'   "1017", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
#'   "1017", "ADJ", NA, "ADVERSE EVENT", "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
#'   "1017", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02")
#' ) %>%
#'   mutate(ASTDTM = ymd_hms(paste(ASTDT, "00:00:00")), AENDTM = ymd_hms(paste(AENDT, "00:00:00")))
#'
#' # Cumulative dose
#' adex %>%
#'   derive_param_exposure(
#'     by_vars = vars(USUBJID),
#'     set_values_to = vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL"),
#'     input_code = "DOSE",
#'     analysis_var = AVAL,
#'     summary_fun = function(x) sum(x, na.rm = TRUE)
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
#'
#' # average dose in w2-24
#' adex %>%
#'   derive_param_exposure(
#'     by_vars = vars(USUBJID),
#'     filter = VISIT %in% c("WEEK 2", "WEEK 24"),
#'     set_values_to = vars(PARAMCD = "AVDW224", PARCAT1 = "WEEK2-24"),
#'     input_code = "DOSE",
#'     analysis_var = AVAL,
#'     summary_fun = function(x) mean(x, na.rm = TRUE)
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
#'
#' # Any dose adjustment?
#' adex %>%
#'   derive_param_exposure(
#'     by_vars = vars(USUBJID),
#'     set_values_to = vars(PARAMCD = "TADJ", PARCAT1 = "OVERALL"),
#'     input_code = "ADJ",
#'     analysis_var = AVALC,
#'     summary_fun = function(x) if_else(sum(!is.na(x)) > 0, "Y", NA_character_)
#'   ) %>%
#'   select(-ASTDTM, -AENDTM)
derive_param_exposure <- function(dataset,
                                  by_vars,
                                  input_code,
                                  analysis_var,
                                  summary_fun,
                                  filter = NULL,
                                  set_values_to = NULL) {
  by_vars <- assert_vars(by_vars)
  analysis_var <- assert_symbol(enquo(analysis_var))

  dtm <- c("ASTDTM", "AENDTM") %in% colnames(dataset)
  dt <- c("ASTDT", "AENDT") %in% colnames(dataset)
  if (all(dtm)) {
    dates <- vars(ASTDTM, AENDTM)
  } else {
    dates <- vars(ASTDT, AENDT)
  }

  assert_data_frame(dataset,
    required_vars = quo_c(by_vars, analysis_var, vars(PARAMCD), dates)
  )
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(input_code)
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(input_code, values = params_available)
  assert_s3_class(summary_fun, "function")

  subset_ds <- dataset %>%
    filter_if(filter)

  add_data <- subset_ds %>%
    get_summary_records(
      by_vars = by_vars,
      filter = PARAMCD == input_code,
      analysis_var = !!analysis_var,
      summary_fun = summary_fun,
      set_values_to = set_values_to
    )

  # add the dates for the derived parameters
  tmp_start <- get_new_tmp_var(dataset)
  tmp_end <- get_new_tmp_var(dataset)
  if (all(dtm)) {
    dates <- subset_ds %>%
      group_by(!!!by_vars) %>%
      summarise(
        !!tmp_start := min(ASTDTM, na.rm = TRUE),
        !!tmp_end := max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
      ) %>%
      ungroup()
    expo_data <- add_data %>%
      derive_vars_merged(dataset_add = dates, by_vars = by_vars) %>%
      mutate(
        ASTDTM = !!tmp_start,
        AENDTM = !!tmp_end
      ) %>%
      remove_tmp_vars()

    if (all(dt)) {
      expo_data <- expo_data %>%
        mutate(ASTDT = date(ASTDTM), AENDT = date(AENDTM))
    }
  } else {
    dates <- subset_ds %>%
      group_by(!!!by_vars) %>%
      summarise(
        !!tmp_start := min(ASTDT, na.rm = TRUE),
        !!tmp_end := max(coalesce(AENDT, ASTDT), na.rm = TRUE)
      ) %>%
      ungroup()
    expo_data <- add_data %>%
      derive_vars_merged(dataset_add = dates, by_vars = by_vars) %>%
      mutate(
        ASTDT = !!tmp_start,
        AENDT = !!tmp_end
      ) %>%
      remove_tmp_vars()
  }

  bind_rows(dataset, expo_data)
}
