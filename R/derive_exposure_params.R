#' Adds a Parameter Computed from the Aggregated Analysis Value of another Parameter
#'
#' Adds a parameter computed from the aggregated analysis value of another parameter.
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
#' @param filter_rows Filter condition
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param new_param Required parameter code
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#' *Permitted Values:* A character value
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
#'   observations. For example `vars(PARCAT1 = "OVERALL")` defines the parameter category
#'   for the new parameter.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset if the
#'   filtered input dataset contains exactly one observation for the parameter
#'   code specified for `input_param`.
#'
#'   For the new observations, `AVAL` is set to the value specified by
#'   `fns` and the variables specified for `set_values_to` are set to
#'   the provided values. The values of the other variables of the input dataset
#'   are retained if they are constant within each by group. Otherwise they are
#'   set to `NA`.
#'
#' @author Samia Kabi
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation bds
#'
#' @export
#'
#' @examples
#'
derive_exposure_params <- function(dataset,
                                   by_vars,
                                   new_param,
                                   input_param,
                                   fns,
                                   filter_rows = NULL,
                                   set_values_to = NULL,
                                   drop_values_from = NULL) {
  assert_data_frame(dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL, AVALC, ASTDTM, AENDTM))
  )
  by_vars <- assert_vars(by_vars)
  assert_character_scalar(new_param)
  assert_character_scalar(input_param)
  filter_rows <- assert_filter_cond(enquo(filter_rows), optional = TRUE)
  assert_vars(drop_values_from, optional = TRUE)
  assert_varval_list(set_values_to, optional = TRUE)
  assert_param_does_not_exist(dataset, new_param)
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(input_param, values = params_available)

  subset_ds <- dataset %>%
    filter_if(filter_rows)

  add_data <- subset_ds %>%
    filter(PARAMCD == input_param) %>%
    derive_summary_records(
      by_vars = by_vars,
      fns = fns,
      set_values_to = vars(PARAMCD___ = !!new_param),
      drop_values_from = drop_values_from
    ) %>%
    filter(PARAMCD___ == new_param) %>%
    mutate(PARAMCD = PARAMCD___) %>%
    # TO UPDATE
    # not sure why the derive_summary_records() fns render a AVAL.x and AVAL.y...
    # when i compute summary for TPDOSE
    # If i remove the summary for TPDOSE it works fine...
    # AVAL.x retain the value of the input param, AVAL.y has the correct result
    select(-ends_with("___"), -ends_with(".x")) %>%
    mutate(
      AVAL = coalesce(!!!select(
        ., starts_with("AVAL"),
        -ends_with("C"),
        -ends_with("C.y")
      )),
      AVALC = coalesce(!!!select(., starts_with("AVAL") &
        (ends_with("C") | ends_with("C.y"))))
    ) %>%
    select(-ends_with(".y"))

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
