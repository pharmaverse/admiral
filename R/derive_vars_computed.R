#' Adds a Variable Computed from the Analysis Value of one or more Parameters
#'
#' Adds a variable computed from the analysis value of one or more parameters.
#' It is expected that the value of the new variable is defined by an expression
#' using the analysis values of other parameters. For example Body Mass Index at
#' Baseline (BMIBL) can be derived from of HEIGHT and WEIGHT parameters.
#'
#' @param dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   The variable specified by `by_vars` must be a unique key of the input
#'   dataset after restricting it by the filter condition (`filter` parameter).
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the additional dataset after restricting it by the filter condition
#'   (`filter_add` parameter) and to the parameters specified by `parameters`.
#'
#'   If the argument is specified, the observations of the additional dataset
#'   are considered in addition to the observations from the input dataset
#'   (`dataset` restricted by `filter`).
#'
#' @param filter Filter condition of dataset
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new variable, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param filter_add Filter condition of additional dataset
#'
#'   The specified condition is applied to the additional dataset before
#'   deriving the new variable, i.e., only observations fulfilling the
#'   condition are taken into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param parameters Required parameter codes
#'
#'   It is expected that all parameter codes (`PARAMCD`) which are required to
#'   derive the new variable are specified for this parameter or the
#'   `constant_parameters` parameter.
#'
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression defines the condition for
#'   selecting the records. For example,
#'   `parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the observations
#'   with `VSTESTCD == "HEIGHT"` from the input data (`dataset` and
#'   `dataset_add`), sets `PARAMCD = "HGHT"` for these observations, and adds
#'   them to the observations to consider.
#'
#'   Unnamed elements in the list of expressions are considered as parameter
#'   codes. For example, `parameters = exprs(WEIGHT, HGHT = VSTESTCD ==
#'   "HEIGHT")` uses the parameter code `"WEIGHT"` and creates a temporary
#'   parameter code `"HGHT"`.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#'
#' @param by_vars Grouping variables
#'
#'   Grouping variables uniquely identifying a set of records for which
#'   `new_vars` are to be calculated.
#'
#'   *Permitted Values:* list of variables
#'
#' @param constant_parameters Required constant parameter codes
#'
#'   It is expected that all the parameter codes (`PARAMCD`) which are required
#'   to derive the new variable and are measured only once are specified here.
#'   For example if BMI should be derived and height is measured only once while
#'   weight is measured at each visit. Height could be specified in the
#'   `constant_parameters` parameter. (Refer to the Example)
#'
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression defines the condition for
#'   selecting the records. For example `constant_parameters =
#'   exprs(HGHT = VSTESTCD == "HEIGHT")` selects the observations with
#'   `VSTESTCD == "HEIGHT"` from the input data (`dataset` and `dataset_add`),
#'   sets `PARAMCD = "HGHT"` for these observations, and adds them to the
#'   observations to consider.
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
#'   to the other parameters using the specified variables.
#'   (Refer to the Example)
#'
#'   *Permitted Values:* list of variables
#'
#'
#' @param new_vars Name of the newly created variables
#'
#'   The specified variables are set to the specified values. The values of
#'   variables of the parameters specified by `parameters` can be accessed using
#'   `<variable name>.<parameter code>`. For example
#'   ```
#'   exprs(
#'     BMIBL = (AVAL.WEIGHT / (AVAL.HEIGHT/100)^2)
#'   )
#'   ```
#'   defines the analysis value and parameter code for the new variable.
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
#'   `by_vars` parameter), the value of the new variable is computed in the
#'   output dataset if the filtered input dataset (`dataset`) or the additional
#'   dataset (`dataset_add`) contains exactly one observation for each parameter
#'   code specified for `parameters`.
#'
#' @return The input dataset with the new variable added.
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
#'
#' # Example 1: Derive BMIBL
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01", "01-1302",   61, "YEARS",
#'   "PILOT01", "17-1344",   64, "YEARS"
#' )
#'
#' advs <- tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
#'   "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", "N",
#'   "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", "N"
#' )
#'
#' derive_vars_computed(
#'   dataset = adsl,
#'   dataset_add = advs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   parameters = c("WEIGHT"),
#'   constant_by_vars = exprs(STUDYID, USUBJID),
#'   constant_parameters = c("HEIGHT"),
#'   new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
#'   filter_add = ABLFL == "Y"
#' )
derive_vars_computed <- function(dataset = NULL,
                                 dataset_add = NULL,
                                 by_vars,
                                 parameters,
                                 new_vars,
                                 filter = NULL,
                                 filter_add = NULL,
                                 constant_by_vars = NULL,
                                 constant_parameters = NULL,
                                 keep_nas = FALSE) {
  assert_vars(by_vars)
  assert_vars(constant_by_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars, optional = TRUE)
  assert_data_frame(dataset_add, optional = TRUE)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)

  assert_varval_list(new_vars)
  if (!is.null(new_vars$PARAMCD) && !is.null(dataset)) {
    assert_param_does_not_exist(dataset, new_vars$PARAMCD)
  }
  assert_logical_scalar(keep_nas)


  parameters <- assert_parameters_argument(parameters)
  constant_parameters <- assert_parameters_argument(constant_parameters, optional = TRUE)

  # select observations and variables required for new observations
  if (is.null(dataset)) {
    data_source <- dataset_add %>%
      filter_if(filter_add)
  } else {
    if (!is.null(dataset_add)) {
      dataset_add <- dataset_add %>%
        filter_if(filter_add)
      data_source <- dataset %>%
        filter_if(filter) %>%
        left_join(dataset_add,
          by = vars2chr(by_vars)
        )
    } else {
      data_source <- dataset %>%
        filter_if(filter)
    }
  }

  temp_return <- get_temp_data(
    data_source,
    by_vars = by_vars,
    parameters = parameters,
    set_values_to = new_vars,
    filter = !!filter,
    filter_add = !!filter_add
  )
  temp_data <- temp_return[["temp_data"]]
  if (is.null(temp_data)) {
    return(dataset)
  }
  analysis_vars_chr <- temp_return[["analysis_vars_chr"]]

  if (!is.null(constant_parameters)) {
    temp_const_data <- get_temp_data(
      data_source,
      by_vars = constant_by_vars,
      parameters = constant_parameters,
      set_values_to = new_vars,
      filter = !!filter,
      filter_add = !!filter_add
    )[["temp_data"]]

    if (is.null(temp_const_data)) {
      return(dataset)
    }

    temp_data <- inner_join(temp_data, temp_const_data, by = vars2chr(constant_by_vars))
  }


  temp_data <- temp_data %>%
    process_set_values_to(new_vars) %>%
    select(-all_of(analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")]))

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  if (!keep_nas) {
    # keep only observations where all analysis values are available
    left_join(dataset,
      temp_data,
      by = vars2chr(by_vars)
    ) %>%
      filter(
        !is.na(!!!parse_exprs(names(new_vars)[1]))
      )
  } else {
    left_join(dataset,
      temp_data,
      by = vars2chr(by_vars)
    )
  }
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
#' @param by_vars By variables
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
#' @param filter_add Filter condition used for restricting the additional dataset
#'
#'    The specified filter condition is used in the warnings only. It is not
#'    applied to the additional dataset.
#'
#'   *Permitted Values:* An unquoted expression
#'
#' @return A dataset with one observation per by group. It contains the
#'   variables specified for `by_vars` and all variables of the form
#'   `<variable>.<parameter>` occurring in `analysis_value`.
#'
#' @keywords internal
get_temp_data <- function(dataset,
                          by_vars,
                          parameters,
                          set_values_to,
                          filter,
                          filter_add) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = by_vars)
  parameters <- assert_parameters_argument(parameters)
  assert_expr_list(set_values_to)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)


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
    if (!is.null(filter) && !is.null(filter_add)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter),
          ") and (",
          expr_label(filter_add),
          ")."
        )
      )
    } else if (is.null(filter)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter_add),
          ")."
        )
      )
    } else if (is.null(filter_add)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter),
          ")."
        )
      )
    }
    return(list(temp_data = NULL))
  }

  params_available <- unique(data_parameters$PARAMCD)
  params_missing <- setdiff(param_values, params_available)
  if (length(params_missing) > 0) {
    if (!is.null(filter) && !is.null(filter_add)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter),
          ") and (",
          expr_label(filter_add),
          ")."
        )
      )
    } else if (is.null(filter)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter_add),
          ")."
        )
      )
    } else if (is.null(filter_add)) {
      warn(
        paste0(
          "The dataset does not contain any observations fullfiling the filter condition (",
          expr_label(filter),
          ")."
        )
      )
    }
    return(list(temp_data = NULL))
  }

  signal_duplicate_records(
    data_parameters,
    by_vars = exprs(!!!by_vars, PARAMCD),
    msg = paste(
      "The filtered dataset contains duplicate records with respect to",
      enumerate(c(vars2chr(by_vars), "PARAMCD")),
      "\nPlease ensure that the variables specified for `by_vars` and `PARAMCD`",
      "are a unique key of the input data set restricted by the condition",
      "specified for `filter` and `filter_add` and to the parameters specified
      for `parameters`."
    )
  )

  # horizontalize data, e.g., AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  analysis_vars <- flatten(map(unname(set_values_to), extract_vars))
  analysis_vars_chr <- vars2chr(analysis_vars)
  multi_dot_names <- str_count(analysis_vars_chr, "\\.") > 1
  if (any(multi_dot_names)) {
    abort(
      paste(
        "The `new_vars` argument contains variable names with more than one dot:",
        enumerate(analysis_vars_chr[multi_dot_names]),
        sep = "\n"
      )
    )
  }
  vars_temp <- analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")] %>%
    str_split(pattern = "\\.") %>%
    map_chr(`[[`, 1) %>%
    unique()

  temp_data <- data_parameters
  for (i in seq_along(vars_temp)) {
    pivoted_data <- pivot_wider(
      select(data_parameters, !!!by_vars, PARAMCD, sym(vars_temp[[i]])),
      names_from = PARAMCD,
      values_from = sym(vars_temp[[i]]),
      names_prefix = paste0(vars_temp[[i]], ".")
    )

    temp_data <- left_join(
      temp_data,
      pivoted_data,
      by = vars2chr(by_vars)
    )
  }

  list(
    temp_data = bind_rows(temp_data) %>%
      select(!!!by_vars, any_of(analysis_vars_chr)),
    analysis_vars_chr = analysis_vars_chr[str_detect(analysis_vars_chr, "\\.")]
  )
}
