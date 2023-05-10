#' Add a Variable Flagging the First or Last Observation Within Each By Group
#'
#' Add a variable flagging the first or last observation within each by group
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param order Sort order
#'
#'   The first or last observation is determined with respect to the specified
#'   order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param new_var Variable to add
#'
#'   The specified variable is added to the output dataset. It is set to `"Y"`
#'   for the first or last observation (depending on the mode) of each by group.
#'
#'   Permitted Values: list of name-value pairs
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
#'   `by_vars` parameter), `new_var` is set to "Y" for the first or last observation
#'   (with respect to the order specified for the `order` parameter and the flag mode
#'   specified for the `mode` parameter). Only observations included by the `filter` parameter
#'   are considered for flagging.
#'   Otherwise, `new_var` is set to `NA`. Thus, the direction of "worst" is considered fixed for
#'   all parameters in the dataset depending on the `order` and the `mode`, i.e. for every
#'   parameter the first or last record will be flagged across the whole dataset.
#'
#' @seealso [derive_var_worst_flag()]
#'
#'
#' @return The input dataset with the new flag variable added
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' example_vs <- tribble(
#'   ~USUBJID, ~VSTESTCD,      ~VISIT, ~VISITNUM, ~VSTPTNUM, ~VSSTRESN,
#'   "1001",     "DIABP", "SCREENING",         1,        10,        64,
#'   "1001",     "DIABP", "SCREENING",         1,        11,        66,
#'   "1001",     "DIABP",  "BASELINE",         2,       100,        68,
#'   "1001",     "DIABP",  "BASELINE",         2,       101,        68,
#'   "1001",     "DIABP",    "WEEK 2",         3,       200,        72,
#'   "1001",     "DIABP",    "WEEK 2",         3,       201,        71,
#'   "1001",     "DIABP",    "WEEK 4",         4,       300,        70,
#'   "1001",     "DIABP",    "WEEK 4",         4,       301,        70
#' )
#'
#' # Flag last value for each patient, test, and visit, baseline observations are ignored
#' example_vs %>%
#'   restrict_derivation(
#'     derivation = derive_var_extreme_flag,
#'     args = params(
#'       by_vars = exprs(USUBJID, VSTESTCD, VISIT),
#'       order = exprs(VSTPTNUM),
#'       new_var = LASTFL,
#'       mode = "last"
#'     ),
#'     filter = VISIT != "BASELINE"
#'   ) %>%
#'   arrange(USUBJID, VSTESTCD, VISITNUM, VSTPTNUM) %>%
#'   select(USUBJID, VSTESTCD, VISIT, VSTPTNUM, VSSTRESN, LASTFL)
#'
#' # Baseline (ABLFL) examples:
#'
#' input <- tribble(
#'   ~STUDYID, ~USUBJID,  ~PARAMCD,     ~AVISIT,                  ~ADT, ~AVAL,    ~DTYPE,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-27"),  15.0,        NA,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-25"),  14.0,        NA,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-23"),  15.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM01",    "WEEK 1", as.Date("2021-04-27"),  10.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM01",    "WEEK 2", as.Date("2021-04-30"),  12.0,        NA,
#'   "TEST01",  "PAT02", "PARAM01", "SCREENING", as.Date("2021-04-27"),  15.0, "AVERAGE",
#'   "TEST01",  "PAT02", "PARAM01",  "BASELINE", as.Date("2021-04-25"),  14.0, "AVERAGE",
#'   "TEST01",  "PAT02", "PARAM01",  "BASELINE", as.Date("2021-04-23"),  15.0, "AVERAGE",
#'   "TEST01",  "PAT02", "PARAM01",    "WEEK 1", as.Date("2021-04-27"),  10.0, "AVERAGE",
#'   "TEST01",  "PAT02", "PARAM01",    "WEEK 2", as.Date("2021-04-30"),  12.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-27"),  15.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-25"),  14.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-23"),  15.0,        NA,
#'   "TEST01",  "PAT01", "PARAM02",  "BASELINE", as.Date("2021-04-27"),  10.0, "AVERAGE",
#'   "TEST01",  "PAT01", "PARAM02",    "WEEK 2", as.Date("2021-04-30"),  12.0,        NA,
#'   "TEST01",  "PAT02", "PARAM02", "SCREENING", as.Date("2021-04-27"),  15.0,        NA,
#'   "TEST01",  "PAT02", "PARAM02",  "BASELINE", as.Date("2021-04-25"),  14.0,        NA,
#'   "TEST01",  "PAT02", "PARAM02",    "WEEK 1", as.Date("2021-04-23"),  15.0,        NA,
#'   "TEST01",  "PAT02", "PARAM02",    "WEEK 1", as.Date("2021-04-27"),  10.0,        NA,
#'   "TEST01",  "PAT02", "PARAM02",  "BASELINE", as.Date("2021-04-30"),  12.0,        NA
#' )
#'
#' # Last observation
#' restrict_derivation(
#'   input,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(ADT),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE"
#' )
#'
#' # Worst observation - Direction = High
#' restrict_derivation(
#'   input,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(AVAL, ADT),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE"
#' )
#'
#' # Worst observation - Direction = Lo
#' restrict_derivation(
#'   input,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(desc(AVAL), ADT),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE"
#' )
#'
#' # Average observation
#' restrict_derivation(
#'   input,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(ADT, desc(AVAL)),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE" & DTYPE == "AVERAGE"
#' )
#'
#' # OCCURDS Examples
#' example_ae <- tribble(
#'   ~USUBJID,         ~AEBODSYS,    ~AEDECOD,   ~AESEV, ~AESTDY, ~AESEQ,
#'   "1015", "GENERAL DISORDERS",  "ERYTHEMA",   "MILD",       2,      1,
#'   "1015", "GENERAL DISORDERS",  "PRURITUS",   "MILD",       2,      2,
#'   "1015",      "GI DISORDERS", "DIARRHOEA",   "MILD",       8,      3,
#'   "1023", "CARDIAC DISORDERS",  "AV BLOCK",   "MILD",      22,      4,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       3,      1,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA", "SEVERE",       5,      2,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       8,      3
#' )
#'
#' # Most severe AE first occurrence per patient
#' example_ae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCIFL,
#'     by_vars = exprs(USUBJID),
#'     order = exprs(TEMP_AESEVN, AESTDY, AESEQ),
#'     mode = "first"
#'   ) %>%
#'   arrange(USUBJID, AESTDY, AESEQ) %>%
#'   select(USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
#'
#' # Most severe AE first occurrence per patient per body system
#' example_ae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCSIFL,
#'     by_vars = exprs(USUBJID, AEBODSYS),
#'     order = exprs(TEMP_AESEVN, AESTDY, AESEQ),
#'     mode = "first"
#'   ) %>%
#'   arrange(USUBJID, AESTDY, AESEQ) %>%
#'   select(USUBJID, AEBODSYS, AESEV, AESTDY, AOCCSIFL)
derive_var_extreme_flag <- function(dataset,
                                    by_vars,
                                    order,
                                    new_var,
                                    mode,
                                    check_type = "warning") {
  new_var <- assert_symbol(enexpr(new_var))
  assert_vars(by_vars)
  assert_expr_list(order)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, !!!extract_vars(order)))
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  check_type <- assert_character_scalar(
    check_type,
    values = c("none", "warning", "error"),
    case_sensitive = FALSE
  )

  # Create flag
  data <- dataset %>%
    derive_var_obs_number(
      new_var = temp_obs_nr,
      order = order,
      by_vars = by_vars,
      check_type = check_type
    )

  if (mode == "first") {
    data <- data %>%
      mutate(!!new_var := if_else(temp_obs_nr == 1, "Y", NA_character_))
  } else {
    data <- data %>%
      group_by(!!!by_vars) %>%
      mutate(!!new_var := if_else(temp_obs_nr == n(), "Y", NA_character_)) %>%
      ungroup()
  }

  # Remove temporary variable
  data %>% select(-temp_obs_nr)
}

#' Adds a Variable Flagging the Maximal / Minimal Value Within a Group of Observations
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*. Please use `slice_derivation()` / `derive_var_extreme_flag()`
#' to derive extreme flags and adjust the `order` argument.
#'
#' @inheritParams derive_var_extreme_flag
#' @param dataset Input dataset.
#' Variables specified by `by_vars`, `order`, `param_var`, and `analysis_var` are expected.
#' @param order Sort order.
#' Used to determine maximal / minimal observation if they are not unique,
#' see Details section for more information.
#' @param new_var Variable to add to the `dataset`.
#' It is set `"Y"` for the maximal / minimal observation of each group,
#' see Details section for more information.
#' @param param_var Variable with the parameter values for which the maximal / minimal
#' value is calculated.
#' @param analysis_var Variable with the measurement values for which the maximal / minimal
#' value is calculated.
#' @param worst_high Character with `param_var` values specifying the parameters
#' referring to "high".
#' Use `character(0)` if not required.
#' @param worst_low Character with `param_var` values specifying the parameters
#' referring to "low".
#' Use `character(0)` if not required.
#'
#' @details For each group with respect to the variables specified by the `by_vars` parameter,
#' the maximal / minimal observation of `analysis_var`
#' is labeled in the `new_var` column as `"Y"`,
#' if its `param_var` is in `worst_high` / `worst_low`.
#' Otherwise, it is assigned `NA`.
#' If there is more than one such maximal / minimal observation,
#' the first one with respect to the order specified by the `order` parameter is flagged. The
#' direction of "worst" depends on the definition of worst for a specified parameters in the
#' arguments `worst_high` / `worst_low`, i.e. for some parameters the highest value is the worst
#' and for others the worst is the lowest value.
#'
#' @seealso [derive_var_extreme_flag()]
#'
#'
#' @return The input dataset with the new flag variable added.
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
derive_var_worst_flag <- function(dataset,
                                  by_vars,
                                  order,
                                  new_var,
                                  param_var,
                                  analysis_var,
                                  worst_high,
                                  worst_low,
                                  check_type = "warning") {
  ### DEPRECATION
  deprecate_stop("0.10.0",
    "derive_var_worst_flag()",
    details = "Please use `slice_derivation()` / `derive_var_extreme_flag()`"
  )
}
