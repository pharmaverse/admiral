#' Adds a Parameter Computed from the Analysis Value of Other Parameters
#'
#' Adds a parameter computed from the analysis value of other parameters. It is
#' expected that the analysis value of the new parameter is defined by an
#' expression using the analysis values of other parameters. For example mean
#' arterial pressure (MAP) can be derived from systolic (SYSBP) and diastolic
#' blood pressure (DIABP) with the formula
#' \deqn{MAP = \frac{SYSBP + 2 DIABP}{3}}{MAP = (SYSBP + 2 DIABP) / 3}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and `AVAL`
#'   are expected.
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
#' @param parameters Required parameter codes
#'
#'   It is expected that all parameter codes (`PARAMCD`) which are required to
#'   derive the new parameter are specified for this parameter or the
#'   `constant_parameters` parameter.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   *Permitted Values:* list of variables
#'
#' @param constant_parameters Required constant parameter codes
#'
#'   It is expected that all the parameter codes (`PARAMCD`) which are required
#'   to derive the new parameter and are measured only once are specified here.
#'   For example if BMI should be derived and height is measured only once while
#'   weight is measured at each visit. Height could be specified in the
#'   `constant_parameters` parameter. (Refer to Example 2)
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values
#'
#' @param constant_by_vars By variables for constant parameters
#'
#'   The constant parameters (parameters that are measured only once) are merged
#'   to the other parameters using the specified variables. (Refer to Example 2)
#'
#'   *Permitted Values:* list of variables
#'
#' @param analysis_value Definition of the analysis value
#'
#'   An expression defining the analysis value (`AVAL`) of the new parameter is
#'   expected. The analysis values of the parameters specified by `parameters`
#'   can be accessed using `AVAL.<parameter code>`, e.g., `AVAL.SYSBP`.
#'
#'   *Permitted Values:* An unquoted expression
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations. For example `vars(PARAMCD = "MAP")` defines the parameter
#'   code for the new parameter.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset if the
#'   filtered input dataset contains exactly one observation for each parameter
#'   code specified for `parameters`.
#'
#'   For the new observations `AVAL` is set to the value specified by
#'   `analysis_value` and the variables specified for `set_values_to` are set to
#'   the provided values. The values of the other variables of the input dataset
#'   are set to `NA`.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' # Example 1: Derive MAP
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
#'   "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
#'   "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
#'   "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
#'   "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
#'   "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
#' )
#'
#' derive_param_computed(
#'   advs,
#'   by_vars = vars(USUBJID, VISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'   set_values_to = vars(
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)",
#'     AVALU = "mmHg"
#'   )
#' )
#'
#' # Example 2: Derive BMI where height is measured only once
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 147, "cm", "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.0, "kg", "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 54.4, "kg", "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)", 53.1, "kg", "WEEK 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 163, "cm", "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 78.5, "kg", "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.3, "kg", "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)", 80.7, "kg", "WEEK 2"
#' )
#'
#' derive_param_computed(
#'   advs,
#'   by_vars = vars(USUBJID, VISIT),
#'   parameters = "WEIGHT",
#'   analysis_value = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
#'   set_values_to = vars(
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)",
#'     AVALU = "kg/m^2"
#'   ),
#'   constant_parameters = c("HEIGHT"),
#'   constant_by_vars = vars(USUBJID)
#' )
derive_param_computed <- function(dataset,
                                  by_vars,
                                  parameters,
                                  analysis_value,
                                  set_values_to,
                                  filter = NULL,
                                  constant_by_vars = NULL,
                                  constant_parameters = NULL) {
  assert_vars(by_vars)
  assert_vars(constant_by_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, PARAMCD, AVAL))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(parameters, values = params_available)
  assert_character_vector(constant_parameters, values = params_available, optional = TRUE)
  assert_varval_list(set_values_to)
  if (!is.null(set_values_to$PARAMCD)) {
    assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  }

  # select observations and variables required for new observations
  data_filtered <- dataset %>%
    filter_if(filter)

  data_parameters <- data_filtered %>%
    filter(PARAMCD %in% parameters)

  if (nrow(data_parameters) == 0L) {
    warn(
      paste0(
        "The input dataset does not contain any observations fullfiling the filter condition (",
        expr_label(filter),
        ") for the parameter codes (PARAMCD) ",
        enumerate(parameters),
        "\nNo new observations were added."
      )
    )
    return(dataset)
  }

  params_available <- unique(data_filtered$PARAMCD)
  params_missing <- setdiff(c(parameters, constant_parameters), params_available)
  if (length(params_missing) > 0) {
    warn(
      paste0(
        "The input dataset does not contain any observations fullfiling the filter condition (",
        expr_label(filter),
        ") for the parameter codes (PARAMCD) ",
        enumerate(params_missing),
        "\nNo new observations were added."
      )
    )
    return(dataset)
  }

  data_parameters <- data_parameters %>%
    select(!!!by_vars, PARAMCD, AVAL)

  signal_duplicate_records(
    data_parameters,
    by_vars = vars(!!!by_vars, PARAMCD),
    msg = paste(
      "The filtered input dataset contains duplicate records with respect to",
      enumerate(c(vars2chr(by_vars), "PARAMCD")),
      "\nPlease ensure that the variables specified for `by_vars` and `PARAMCD`",
      "are a unique key of the input data set restricted by the condition",
      "specified for `filter` and to the parameters specified for `parameters`."
    )
  )

  # horizontalize data, AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  hori_data <- data_parameters %>%
    pivot_wider(names_from = PARAMCD, values_from = AVAL, names_prefix = "AVAL.")

  if (!is.null(constant_parameters)) {
    data_const_parameters <- data_filtered %>%
      filter(PARAMCD %in% constant_parameters) %>%
      select(!!!vars(!!!constant_by_vars, PARAMCD, AVAL))

    hori_const_data <- data_const_parameters %>%
      pivot_wider(names_from = PARAMCD, values_from = AVAL, names_prefix = "AVAL.")

    hori_data <- inner_join(hori_data, hori_const_data, by = vars2chr(constant_by_vars))
  }

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  hori_data <- hori_data %>%
    # keep only observations where all analysis values are available
    filter(!!!parse_exprs(map_chr(
      c(parameters, constant_parameters),
      ~ str_c("!is.na(AVAL.", .x, ")")
    ))) %>%
    mutate(AVAL = !!enquo(analysis_value), !!!set_values_to) %>%
    select(-starts_with("AVAL."))

  bind_rows(dataset, hori_data)
}

#' Adds a Parameter Computed from the Analysis Value of Other Parameters
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated.  Please use `derive_param-computed()` instead.
#'
#' @inheritParams derive_param_computed
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @keywords deprecated
#' @family deprecated
#'
#' @export
#'
derive_derived_param <- function(dataset,
                                 by_vars,
                                 parameters,
                                 analysis_value,
                                 set_values_to,
                                 filter = NULL,
                                 constant_by_vars = NULL,
                                 constant_parameters = NULL) {
  deprecate_warn("0.8.0", "derive_derived_param()", "derive_param_computed()")
  derive_param_computed(
    dataset,
    by_vars = by_vars,
    parameters = parameters,
    analysis_value = !!enquo(analysis_value),
    set_values_to = set_values_to,
    filter = !!enquo(filter),
    constant_by_vars = constant_by_vars,
    constant_parameters = constant_parameters
  )
}
