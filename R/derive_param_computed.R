#' Adds a Parameter Computed from the Analysis Value of Other Parameters
#'
#' Adds a parameter computed from the analysis value of other parameters. It is
#' expected that the analysis value of the new parameter is defined by an
#' expression using the analysis values of other parameters, such as addition/sum,
#' subtraction/difference, multiplication/product, division/ratio,
#' exponentiation/logarithm, or by formula.
#' <br/><br/>
#' For example mean arterial pressure (MAP) can be derived from systolic (SYSBP)
#' and diastolic blood pressure (DIABP) with the formula
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
#' @permitted [dataset]
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
#' @permitted [dataset]
#'
#' @param filter Filter condition
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#' @permitted [condition]
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
#' @permitted A character vector of `PARAMCD` values or a list of expressions
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
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
#' @permitted A character vector of `PARAMCD` values or a list of expressions
#'
#' @param constant_by_vars By variables for constant parameters
#'
#'   The constant parameters (parameters that are measured only once) are merged
#'   to the other parameters using the specified variables. (Refer to Example 2)
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
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
#'   Note that `dplyr` helper functions such as `dplyr::starts_with()` should
#'   be avoided unless the list of variable-value pairs is clearly
#'   specified in a statement via the `set_values_to` argument.
#'
#' @permitted [expr_list_formula]
#'
#' @param keep_nas Keep observations with `NA`s
#'
#'   If the argument is set to `TRUE`, observations are added even if some of
#'   the values contributing to the computed value are `NA` (see Example 1b).
#'
#'   If the argument is set to a list of variables, observations are added even
#'   if some of specified variables are `NA` (see Example 1c).
#'
#' @permitted `TRUE`, `FALSE`, or a list of variables created by
#'   `exprs()` e.g. `exprs(ADTF, ATMF)`
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset if the
#'   filtered input dataset (`dataset`) or the additional dataset
#'   (`dataset_add`) contains exactly one observation for each parameter code
#'   specified for `parameters` and all contributing values like `AVAL.SYSBP`
#'   are not `NA`. The `keep_nas` can be used to specify variables for which
#'   `NA`s are acceptable. See also Example 1b and 1c.
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
#' @examplesx
#'
#' @caption Example 1 - Data setup
#'
#' @info Examples 1a, 1b, and 1c use the following `ADVS` data.
#'
#' @code
#' ADVS <- tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
#'   "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
#'   "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
#'   "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
#'   "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
#' ) %>%
#'   mutate(
#'     AVALU = "mmHg",
#'     ADT = case_when(
#'       VISIT == "BASELINE" ~ as.Date("2024-01-10"),
#'       VISIT == "WEEK 2" ~ as.Date("2024-01-24")
#'     ),
#'     ADTF = NA_character_
#'   )
#'
#' @caption Example 1a - Adding a parameter computed from a formula
#'   (`parameters`, `set_values_to`)
#'
#' @info Derive mean arterial pressure (MAP) from systolic (SYSBP)
#'   and diastolic blood pressure (DIABP).
#'
#' - Here, for each `USUBJID` and `VISIT` group (specified in `by_vars`),
#'   an observation is added to the output dataset when the filtered
#'   input dataset (`dataset`) contains exactly one observation for
#'   each parameter code specified for `parameters` and all contributing
#'   values (e.g., `AVAL.SYSBP` and `AVAL.DIABP`) are not `NA`.
#'   Indeed, patient `01-701-1028` does not get a `"WEEK 2"`-derived record
#'   as `AVAL` is `NA` for their `"WEEK 2"` systolic blood pressure.
#'
#' @code
#' derive_param_computed(
#'   ADVS,
#'   by_vars = exprs(USUBJID, VISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)",
#'     AVALU = "mmHg",
#'     ADT = ADT.SYSBP
#'   )
#' ) %>%
#' select(-PARAM)
#'
#' @caption Example 1b - Keeping missing values for any source
#'     variables (`keep_nas = TRUE`)
#'
#' @info Use option `keep_nas = TRUE` to derive MAP in the case where
#'     some/all values of a variable used in the computation are missing.
#'
#' - Note that observations will be added here even if some of the values contributing
#'   to the computed values are `NA`. In particular, patient `01-701-1028`
#'   does get a `"WEEK 2"`-derived record as compared to Example 1a, but
#'   with `AVAL = NA`.
#'
#' @code
#' derive_param_computed(
#'   ADVS,
#'   by_vars = exprs(USUBJID, VISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)",
#'     AVALU = "mmHg",
#'     ADT = ADT.SYSBP,
#'     ADTF = ADTF.SYSBP
#'   ),
#'   keep_nas = TRUE
#' )%>%
#' select(-PARAM)
#'
#' @caption Example 1c - Keeping missing values for some source
#'     variables (`keep_nas = exprs()`)
#'
#' @info Use option `keep_nas = exprs(ADTF)` to derive MAP in the case where
#'     some/all values of a variable used in the computation are
#'     missing but keeping `NA` values of `ADTF`.
#'
#' - This is subtly distinct from Examples 1a and 1b. In 1a, we do not
#'   get new derived records if any of the source records have a value
#'   of `NA` for a variable that is included in `set_values_to`.
#'   In 1b, we do the opposite and allow the creation of new records
#'   regardless of how many `NA`s we encounter in the source variables.
#' - Here, we want to disregard `NA` values but only from the variables
#'   that are specified via `keep_na_values`.
#' - This is important because we have added `ADTF` in `set_values_to`,
#'   but all values of this variable are `NA`. As such, in order to
#'   get any derived records at all, but continue not getting one
#'   when `AVAL` is `NA` in any of the source records,
#'   (see patient `"01-701-1028"` again), we specify `keep_nas = exprs(ADTF)`.
#'
#' @code
#' derive_param_computed(
#'   ADVS,
#'   by_vars = exprs(USUBJID, VISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP",
#'     PARAM = "Mean Arterial Pressure (mmHg)",
#'     AVALU = "mmHg",
#'     ADT = ADT.SYSBP,
#'     ADTF = ADTF.SYSBP
#'   ),
#'   keep_nas = exprs(ADTF)
#' )
#'
#' @caption Example 2 - Derivations using parameters measured only once
#' (`constant_parameters` and `constant_by_vars`)
#'
#' @info Derive BMI where `HEIGHT` is measured only once.
#'
#' - In the above examples, for each parameter specified in the
#'   `parameters` argument, we expect one record per by group, where the by
#'   group is specified in `by_vars`. However, if a parameter is only
#'   measured once, it can be specified in `constant_parameters` instead.
#' - A modified by group still needs to be provided for the constant
#'   parameters. This can be done via `constant_by_vars`.
#' - See the example below, where weight is measured for each patient
#'   at each visit (`by_vars = exprs(USUBJID, VISIT)`), while height
#'   is measured for each patient only at the first visit
#'   (`constant_parameters = "HEIGHT"`, `constant_by_vars = exprs(USUBJID`)).
#'
#' @code
#' ADVS <- tribble(
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
#'   ADVS,
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
#' @caption Example 3 - Derivations including data from an additional
#' dataset (`dataset_add`) and non-`AVAL` variables
#'
#' @info Use data from an additional dataset and other variables than `AVAL`.
#'
#' - In this example, the dataset specified via `dataset_add` (e.g., `QS`)
#'   is an SDTM dataset. There is no parameter code in the dataset.
#' - The `parameters` argument is therefore used to specify a list of
#'   expressions to derive temporary parameter codes.
#' - Then, `set_values_to` is used to specify the values for the new
#'   observations of each variable, and variable-value pairs from both
#'   datasets are referenced via `exprs()`.
#'
#' @code
#' QS <- tribble(
#'   ~USUBJID, ~AVISIT,   ~QSTESTCD, ~QSORRES, ~QSSTRESN,
#'   "1",      "WEEK 2",  "CHSF112", NA,               1,
#'   "1",      "WEEK 2",  "CHSF113", "Yes",           NA,
#'   "1",      "WEEK 2",  "CHSF114", NA,               1,
#'   "1",      "WEEK 4",  "CHSF112", NA,               2,
#'   "1",      "WEEK 4",  "CHSF113", "No",            NA,
#'   "1",      "WEEK 4",  "CHSF114", NA,               1
#' )
#'
#' ADCHSF <- tribble(
#'   ~USUBJID, ~AVISIT,  ~PARAMCD, ~QSSTRESN, ~AVAL,
#'   "1",      "WEEK 2", "CHSF12", 1,             6,
#'   "1",      "WEEK 2", "CHSF14", 1,             6,
#'   "1",      "WEEK 4", "CHSF12", 2,            12,
#'   "1",      "WEEK 4", "CHSF14", 1,             6
#' ) %>%
#'   mutate(QSORRES = NA_character_)
#'
#' derive_param_computed(
#'   ADCHSF,
#'   dataset_add = QS,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   parameters = exprs(CHSF12, CHSF13 = QSTESTCD %in% c("CHSF113"), CHSF14),
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
#' @caption Example 4 - Computing more than one variable
#'
#' @info Specify more than one variable-value pair via `set_values_to`.
#'
#' - In this example, the values of `AVALC`, `ADTM`, `ADTF`, `PARAMCD`,
#'   and `PARAM` are determined via distinctly defined analysis values
#'   and parameter codes.
#' - This is different from Example 3 as more than one variable is
#'   derived.
#'
#' @code
#' ADLB_TBILIALK <- tribble(
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
#'   dataset_add = ADLB_TBILIALK,
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
  if (typeof(keep_nas) == "list") {
    assert_vars(keep_nas)
  } else {
    assert_logical_scalar(
      keep_nas,
      message = paste(
        "Argument {.arg {arg_name}} must be either {.val {TRUE}}, {.val {FALSE}},",
        "or a list of {.cls symbol}, but is {.obj_type_friendly {arg}}."
      )
    )
  }

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

  if (isFALSE(keep_nas) || typeof(keep_nas) == "list") {
    # keep only observations where the specified analysis values are available
    if (typeof(keep_nas) == "list") {
      na_vars <- discard(
        analysis_vars_chr,
        ~ str_detect(., paste0("^(", paste(vars2chr(keep_nas), collapse = "|"), ")\\."))
      )
    } else {
      na_vars <- analysis_vars_chr
    }
    nobs_before <- nrow(hori_data)
    hori_data <- filter(
      hori_data,
      !!!parse_exprs(map_chr(
        na_vars,
        ~ str_c("!is.na(", .x, ")")
      ))
    )
    if (nobs_before > 0 && nrow(hori_data) == 0) {
      cli_inform(
        c(
          paste(
            "No computed records were added because for all potential computed",
            "records at least one of the contributing values was {.val {NA}}."
          ),
          paste(
            "If this is not expected, please check the input data and the value of",
            "the {.arg keep_nas} argument."
          )
        ),
        class = "derive_param_computed_all_na"
      )
    }
  }

  # add computed variables like AVAL and constant variables like PARAMCD
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
      cli_abort(
        paste(
          "{.arg {rlang::caller_arg(parameters)}} must be a character vector",
          "or a list of expressions but it is {.obj_type_friendly {parameters}}."
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
#' @permitted A character vector of `PARAMCD` values or a list of expressions
#'
#' @param set_values_to
#'
#'   All variables of the form `<variable>.<parameter>` like `AVAL.WEIGHT` are
#'   added to the input dataset. They are set to the value of the variable for
#'   the parameter. E.g., `AVAL.WEIGHT` is set to the value of `AVAL` where
#'   `PARAMCD == "WEIGHT"`.
#'
#' @permitted
#'
#' @param filter Filter condition used for restricting the input dataset
#'
#'    The specified filter condition is used in the warnings only. It is not
#'    applied to the input dataset.
#'
#' @permitted An unquoted expression
#'
#' @return A dataset with one observation per by group. It contains the
#'   variables specified for `by_vars` and all variables of the form
#'   `<variable>.<parameter>` occurring in `set_values_to`.
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
    cli_warn(
      c(paste0(
        "The input dataset does not contain any observations fullfiling the filter condition (",
        "{.code {expr_label(filter)}}}",
        ") for the parameter codes (PARAMCD) ",
        "{.val {param_values}}",
        i = "No new observations were added."
      ))
    )
    return(list(hori_data = NULL))
  }

  params_available <- unique(data_parameters$PARAMCD)
  params_missing <- setdiff(param_values, params_available)
  if (length(params_missing) > 0) {
    cli_warn(
      paste0(
        "The input dataset does not contain any observations fullfiling the filter condition (",
        "{.code {expr_label(filter)}}",
        ") for the parameter codes (PARAMCD) ",
        "{.val {params_missing}}",
        i = "No new observations were added."
      )
    )
    return(list(hori_data = NULL))
  }

  signal_duplicate_records(
    data_parameters,
    by_vars = exprs(!!!by_vars, PARAMCD),
    msg = c(
      paste(
        "The filtered input dataset contains duplicate records with respect to",
        "{.var {c(vars2chr(by_vars), \"PARAMCD\")}}"
      ),
      i = "Please ensure that the variables specified for {.arg by_vars} and {.var PARAMCD}",
      "are a unique key of the input data set restricted by the condition",
      "specified for {.arg filter} and to the parameters specified for {.arg parameters}."
    )
  )

  # horizontalize data, e.g., AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  analysis_vars <- flatten(map(unname(set_values_to), extract_vars))
  analysis_vars_chr <- vars2chr(analysis_vars)
  multi_dot_names <- str_count(analysis_vars_chr, "\\.") > 1
  if (any(multi_dot_names)) {
    cli_abort(
      c(
        "The `set_values_to` argument contains variable names with more than one dot:",
        "{.var {analysis_vars_chr[multi_dot_names]}}"
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
