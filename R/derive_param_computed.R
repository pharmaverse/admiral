#' Adds a Parameter Computed from the Analysis Value of Other Parameters
#'
#' Adds a parameter computed from the analysis value of other parameters. It is
#' expected that the analysis value of the new parameter is defined by an
#' expression using the analysis values of other parameters. For example mean
#' arterial pressure (MAP) can be derived from systolic (SYSBP) and diastolic
#' blood pressure (DIABP) with the formula
#' \deqn{MAP = \frac{SYSBP + 2 DIABP}{3}}{MAP = (SYSBP + 2 DIABP) / 3}
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'  `PARAMCD` is expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `parameters`.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the additional dataset after restricting it to the parameters specified by
#'   `parameters`.
#'
#'   If the argument is specified, the observations of the additional dataset
#'   are considered in addition to the observations from the input dataset
#'   (`dataset` restricted by `filter`).
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
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression the condition for selecting the
#'   records. For example `parameters = exprs(HGHT = VSTESTCD == "HEIGHT")`
#'   selects the observations with `VSTESTCD == "HEIGHT"` from the input data
#'   (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
#'   observations, and adds them to the observations to consider.
#'
#'   Unnamed elements in the list of expressions are considered as parameter
#'   codes. For example, `parameters = exprs(WEIGHT, HGHT = VSTESTCD ==
#'   "HEIGHT")` uses the parameter code `"WEIGHT"` and creates a temporary
#'   parameter code `"HGHT"`.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#' @param analysis_var Analysis variable
#'
#'   `r lifecycle::badge("deprecated")` Please use `set_values_to` instead.
#'
#'   The specified variable is set to the value of `analysis_value` for the new
#'   observations.
#'
#'   *Permitted Values*: An unquoted symbol
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param constant_parameters Required constant parameter codes
#'
#'   It is expected that all the parameter codes (`PARAMCD`) which are required
#'   to derive the new parameter and are measured only once are specified here.
#'   For example if BMI should be derived and height is measured only once while
#'   weight is measured at each visit. Height could be specified in the
#'   `constant_parameters` parameter. (Refer to Example 2)
#'
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression the condition for selecting the
#'   records. For example `constant_parameters = exprs(HGHT = VSTESTCD ==
#'   "HEIGHT")` selects the observations with `VSTESTCD == "HEIGHT"` from the
#'   input data (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
#'   observations, and adds them to the observations to consider.
#'
#'   Unnamed elements in the list of expressions are considered as parameter
#'   codes. For example, `constant_parameters = exprs(WEIGHT, HGHT = VSTESTCD ==
#'   "HEIGHT")` uses the parameter code `"WEIGHT"` and creates a temporary
#'   parameter code `"HGHT"`.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#' @param constant_by_vars By variables for constant parameters
#'
#'   The constant parameters (parameters that are measured only once) are merged
#'   to the other parameters using the specified variables. (Refer to Example 2)
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param analysis_value Definition of the analysis value
#'
#'  `r lifecycle::badge("deprecated")` Please use `set_values_to` instead.
#'
#'   An expression defining the analysis value (`AVAL`) of the new parameter is
#'   expected. The values of variables of the parameters specified by
#'   `parameters` can be accessed using `<variable name>.<parameter code>`,
#'   e.g., `AVAL.SYSBP`.
#'
#'   Variable names in the expression must not contain more than one dot.
#'
#'   *Permitted Values:* An unquoted expression
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations. The values of variables of the parameters specified by
#'   `parameters` can be accessed using `<variable name>.<parameter code>`. For
#'   example
#'   ```
#'   exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP"
#'   )
#'   ```
#'   defines the analysis value and parameter code for the new parameter.
#'
#'   Variable names in the expression must not contain more than one dot.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @param keep_nas Keep observations with `NA`s
#'
#'   If the argument is set to `TRUE`, observations are added even if some of
#'   the values contributing to the computed value are `NA`.
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset if the
#'   filtered input dataset (`dataset`) or the additional dataset
#'   (`dataset_add`) contains exactly one observation for each parameter code
#'   specified for `parameters`.
#'
#'   For the new observations the variables specified for `set_values_to` are
#'   set to the provided values. The values of the other variables of the input
#'   dataset are set to `NA`.
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
#' library(dplyr)
#' library(lubridate)
#'
#' # Example 1: Derive MAP
#' advs <- tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "mmHg", "BASELINE",
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "mmHg", "WEEK 2",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "mmHg", "BASELINE",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "mmHg", "WEEK 2",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "mmHg", "BASELINE",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "mmHg", "WEEK 2",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "mmHg", "BASELINE",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    132, "mmHg", "WEEK 2"
#' )
#'
#' derive_param_computed(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)",
#'     AVALU = "mmHg"
#'   )
#' )
#'
#' # Example 2: Derive BMI where height is measured only once
#' advs <- tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HEIGHT", "Height (cm)", 147.0, "cm",   "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
#'   "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
#'   "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
#'   "01-701-1028", "HEIGHT", "Height (cm)", 163.0, "cm",   "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
#'   "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
#'   "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
#' )
#'
#' derive_param_computed(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   parameters = "WEIGHT",
#'   set_values_to = exprs(
#'     AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
#'     PARAMCD = "BMI",
#'     PARAM = "Body Mass Index (kg/m^2)",
#'     AVALU = "kg/m^2"
#'   ),
#'   constant_parameters = c("HEIGHT"),
#'   constant_by_vars = exprs(USUBJID)
#' )
#'
#' # Example 3: Using data from an additional dataset and other variables than AVAL
#' qs <- tribble(
#'   ~USUBJID, ~AVISIT,   ~QSTESTCD, ~QSORRES, ~QSSTRESN,
#'   "1",      "WEEK 2",  "CHSF112", NA,               1,
#'   "1",      "WEEK 2",  "CHSF113", "Yes",           NA,
#'   "1",      "WEEK 2",  "CHSF114", NA,               1,
#'   "1",      "WEEK 4",  "CHSF112", NA,               2,
#'   "1",      "WEEK 4",  "CHSF113", "No",            NA,
#'   "1",      "WEEK 4",  "CHSF114", NA,               1
#' )
#'
#' adchsf <- tribble(
#'   ~USUBJID, ~AVISIT,  ~PARAMCD, ~QSSTRESN, ~AVAL,
#'   "1",      "WEEK 2", "CHSF12", 1,             6,
#'   "1",      "WEEK 2", "CHSF14", 1,             6,
#'   "1",      "WEEK 4", "CHSF12", 2,            12,
#'   "1",      "WEEK 4", "CHSF14", 1,             6
#' ) %>%
#'   mutate(QSORRES = NA_character_)
#'
#' derive_param_computed(
#'   adchsf,
#'   dataset_add = qs,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   parameters = exprs(CHSF12, CHSF13 = QSTESTCD %in% c("CHSF113", "CHSF213"), CHSF14),
#'   set_values_to = exprs(
#'     AVAL = case_when(
#'       QSORRES.CHSF13 == "Not applicable" ~ 0,
#'       QSORRES.CHSF13 == "Yes" ~ 38,
#'       QSORRES.CHSF13 == "No" ~ if_else(
#'         QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
#'         25,
#'         0
#'       )
#'     ),
#'     PARAMCD = "CHSF13"
#'   )
#' )
#'
#' # Example 4: Computing more than one variable
#' adlb_tbilialk <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVALC, ~ADTM,        ~ADTF,
#'   "1",      "ALK2",   "Y",    "2021-05-13", NA_character_,
#'   "1",      "TBILI2", "Y",    "2021-06-30", "D",
#'   "2",      "ALK2",   "Y",    "2021-12-31", "M",
#'   "2",      "TBILI2", "N",    "2021-11-11", NA_character_,
#'   "3",      "ALK2",   "N",    "2021-04-03", NA_character_,
#'   "3",      "TBILI2", "N",    "2021-04-04", NA_character_
#' ) %>%
#'   mutate(ADTM = ymd(ADTM))
#'
#' derive_param_computed(
#'   dataset_add = adlb_tbilialk,
#'   by_vars = exprs(USUBJID),
#'   parameters = c("ALK2", "TBILI2"),
#'   set_values_to = exprs(
#'     AVALC = if_else(AVALC.TBILI2 == "Y" & AVALC.ALK2 == "Y", "Y", "N"),
#'     ADTM = pmax(ADTM.TBILI2, ADTM.ALK2),
#'     ADTF = if_else(ADTM == ADTM.TBILI2, ADTF.TBILI2, ADTF.ALK2),
#'     PARAMCD = "TB2AK2",
#'     PARAM = "TBILI > 2 times ULN and ALKPH <= 2 times ULN"
#'   ),
#'   keep_nas = TRUE
#' )
derive_param_computed <- function(dataset = NULL,
                                  dataset_add = NULL,
                                  by_vars,
                                  parameters,
                                  analysis_var = AVAL,
                                  analysis_value,
                                  set_values_to,
                                  filter = NULL,
                                  constant_by_vars = NULL,
                                  constant_parameters = NULL,
                                  keep_nas = FALSE) {
  assert_vars(by_vars)
  assert_vars(constant_by_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars, optional = TRUE)
  assert_data_frame(dataset_add, optional = TRUE)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_varval_list(set_values_to)
  if (!is.null(set_values_to$PARAMCD) && !is.null(dataset)) {
    assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  }
  assert_logical_scalar(keep_nas)
  ### BEGIN DEPRECATION
  if (!missing(analysis_var)) {
    deprecate_stop(
      "0.12.0",
      "derive_param_computed(analysis_var = )",
      "derive_param_computed(set_values_to = )"
    )
  }
  analysis_var <- assert_symbol(enexpr(analysis_var))

  if (!missing(analysis_value)) {
    deprecate_stop(
      "0.12.0",
      "derive_param_computed(analysis_value = )",
      "derive_param_computed(set_values_to = )"
    )
    set_values_to <- exprs(!!analysis_var := !!enexpr(analysis_value), !!!set_values_to)
  }
  ### END DEPRECATION

  parameters <- assert_parameters_argument(parameters)
  constant_parameters <- assert_parameters_argument(constant_parameters, optional = TRUE)

  # select observations and variables required for new observations
  if (is.null(dataset)) {
    data_source <- dataset_add
  } else {
    data_source <- dataset %>%
      filter_if(filter) %>%
      bind_rows(dataset_add)
  }

  hori_return <- get_hori_data(
    data_source,
    by_vars = by_vars,
    parameters = parameters,
    set_values_to = set_values_to,
    filter = !!filter
  )
  hori_data <- hori_return[["hori_data"]]
  if (is.null(hori_data)) {
    return(dataset)
  }
  analysis_vars_chr <- hori_return[["analysis_vars_chr"]]

  if (!is.null(constant_parameters)) {
    hori_const_data <- get_hori_data(
      data_source,
      by_vars = constant_by_vars,
      parameters = constant_parameters,
      set_values_to = set_values_to,
      filter = !!filter
    )[["hori_data"]]

    if (is.null(hori_const_data)) {
      return(dataset)
    }

    hori_data <- inner_join(hori_data, hori_const_data, by = vars2chr(constant_by_vars))
  }

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  if (!keep_nas) {
    # keep only observations where all analysis values are available
    hori_data <- filter(
      hori_data,
      !!!parse_exprs(map_chr(
        analysis_vars_chr,
        ~ str_c("!is.na(", .x, ")")
      ))
    )
  }
  hori_data <- hori_data %>%
    process_set_values_to(set_values_to) %>%
    select(-all_of(analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")]))

  bind_rows(dataset, hori_data)
}

#' Asserts `parameters` Argument and Converts to List of Expressions
#'
#' The function asserts that the argument is a character vector or a list of
#' expressions. If it is a character vector, it converts it to a list of
#' symbols.
#'
#' @param parameters The argument to check
#'
#' @param optional Is the checked argument optional? If set to `FALSE` and
#'   `parameters` is `NULL` then an error is thrown.
#'
#' @return The `parameters` argument (converted to a list of symbol, if it is a
#'   character vector)
#'
#' @keywords internal
assert_parameters_argument <- function(parameters, optional = TRUE) {
  assert_logical_scalar(optional)
  if (optional && is.null(parameters)) {
    return(invisible(parameters))
  }

  if (typeof(parameters) == "character") {
    parameters <- map(parameters, sym)
  } else {
    if (!inherits(parameters, "list") || any(!map_lgl(
      parameters,
      ~ is_call(.x) || is_expression(.x)
    ))) {
      abort(
        paste0(
          "`",
          arg_name(substitute(parameters)),
          "` must be a character vector or a list of expressions but it is ",
          what_is_it(parameters),
          "."
        )
      )
    }
  }
  parameters
}

#' Creating Temporary Parameters and `<variable>.<parameter>` Variables
#'
#' The function creates temporary parameters and variables of the form
#' `<variable>.<parameter>`, e.g., `AVAL.WEIGHT`.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param by_vars Grouping variables
#'
#' `r roxygen_param_by_vars()`
#'
#' @param parameters List of parameter codes
#'
#'   The input dataset is restricted to the specified parameter codes. If an
#'   expression is specified, a new parameter code is added to the input
#'   dataset. The name of the element defines the parameter code and the
#'   expression the observations to select.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#' @param set_values_to
#'
#'   All variables of the form `<variable>.<parameter>` like `AVAL.WEIGHT` are
#'   added to the input dataset. They are set to the value of the variable for
#'   the parameter. E.g., `AVAL.WEIGHT` is set to the value of `AVAL` where
#'   `PARAMCD == "WEIGHT"`.
#'
#'   *Permitted Values:* A list of expressions
#'
#' @param filter Filter condition used for restricting the input dataset
#'
#'    The specified filter condition is used in the warnings only. It is not
#'    applied to the input dataset.
#'
#'   *Permitted Values:* An unquoted expression
#'
#' @return A dataset with one observation per by group. It contains the
#'   variables specified for `by_vars` and all variables of the form
#'   `<variable>.<parameter>` occurring in `analysis_value`.
#'
#' @keywords internal
get_hori_data <- function(dataset,
                          by_vars,
                          parameters,
                          set_values_to,
                          filter) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = by_vars)
  parameters <- assert_parameters_argument(parameters)
  assert_expr_list(set_values_to)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  # determine parameter values
  if (is.null(names(parameters))) {
    param_values <- map(parameters, as_label)
  } else {
    param_values <- map2(
      parameters,
      names(parameters),
      ~ if_else(.y == "", as_label(.x), .y)
    )
  }

  new_params <- parameters[names(parameters) != ""]
  new_names <- names(new_params)

  new_data <- vector("list", length(new_params))
  for (i in seq_along(new_params)) {
    new_data[[i]] <- filter(dataset, !!new_params[[i]]) %>%
      mutate(PARAMCD = new_names[[i]])
  }

  data_parameters <- dataset %>%
    bind_rows(new_data) %>%
    filter(PARAMCD %in% param_values)

  if (nrow(data_parameters) == 0L) {
    warn(
      paste0(
        "The input dataset does not contain any observations fullfiling the filter condition (",
        expr_label(filter),
        ") for the parameter codes (PARAMCD) ",
        enumerate(param_values),
        "\nNo new observations were added."
      )
    )
    return(list(hori_data = NULL))
  }

  params_available <- unique(data_parameters$PARAMCD)
  params_missing <- setdiff(param_values, params_available)
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
    return(list(hori_data = NULL))
  }

  signal_duplicate_records(
    data_parameters,
    by_vars = exprs(!!!by_vars, PARAMCD),
    msg = paste(
      "The filtered input dataset contains duplicate records with respect to",
      enumerate(c(vars2chr(by_vars), "PARAMCD")),
      "\nPlease ensure that the variables specified for `by_vars` and `PARAMCD`",
      "are a unique key of the input data set restricted by the condition",
      "specified for `filter` and to the parameters specified for `parameters`."
    )
  )

  # horizontalize data, e.g., AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  analysis_vars <- flatten(map(unname(set_values_to), extract_vars))
  analysis_vars_chr <- vars2chr(analysis_vars)
  multi_dot_names <- str_count(analysis_vars_chr, "\\.") > 1
  if (any(multi_dot_names)) {
    abort(
      paste(
        "The `set_values_to` argument contains variable names with more than on dot:",
        enumerate(analysis_vars_chr[multi_dot_names]),
        sep = "\n"
      )
    )
  }
  vars_hori <- analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")] %>%
    str_split(pattern = "\\.") %>%
    map_chr(`[[`, 1) %>%
    unique()

  hori_data <- data_parameters
  for (i in seq_along(vars_hori)) {
    pivoted_data <- pivot_wider(
      select(data_parameters, !!!by_vars, PARAMCD, sym(vars_hori[[i]])),
      names_from = PARAMCD,
      values_from = sym(vars_hori[[i]]),
      names_prefix = paste0(vars_hori[[i]], ".")
    )
    if (i == 1) {
      hori_data <- pivoted_data
    } else {
      hori_data <- left_join(
        hori_data,
        pivoted_data,
        by = vars2chr(by_vars)
      )
    }
  }

  list(
    hori_data = bind_rows(hori_data) %>%
      select(!!!by_vars, any_of(analysis_vars_chr)),
    analysis_vars_chr = analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")]
  )
}
