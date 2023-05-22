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
#' library(dplyr, warn.conflicts = FALSE)
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~VSTESTCD, ~VSSTRESN,      ~VISIT, ~VSTPTNUM,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        70, "SCREENING",       815,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        78, "SCREENING",       816,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        80, "SCREENING",       817,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        68,  "BASELINE",       815,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        78,  "BASELINE",       816,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        74,  "BASELINE",       817,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        60,    "WEEK 4",       815,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        70,    "WEEK 4",       816,
#'   "PILOT01",    "VS", "10-1187",   "DIABP",        74,    "WEEK 4",       817,
#'   "PILOT01",    "VS", "10-1187",  "WEIGHT",      49.9, "SCREENING",        NA,
#'   "PILOT01",    "VS", "10-1187",  "WEIGHT",      49.9,  "BASELINE",        NA,
#'   "PILOT01",    "VS", "10-1187",  "WEIGHT",     49.44,    "WEEK 4",        NA,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        68, "SCREENING",       815,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        70, "SCREENING",       816,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        72, "SCREENING",       817,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        80,  "BASELINE",       815,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        86,  "BASELINE",       816,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        84,  "BASELINE",       817,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        80,    "WEEK 4",       815,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        88,    "WEEK 4",       816,
#'   "PILOT01",    "VS", "10-1271",   "DIABP",        86,    "WEEK 4",       817,
#'   "PILOT01",    "VS", "10-1271",  "WEIGHT",     47.17, "SCREENING",        NA,
#'   "PILOT01",    "VS", "10-1271",  "WEIGHT",     47.63,  "BASELINE",        NA,
#'   "PILOT01",    "VS", "10-1271",  "WEIGHT",     49.44,    "WEEK 4",        NA
#' )
#'
#' # Flag last value for each patient, test, and visit, baseline observations are ignored
#' vs %>%
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
#'   arrange(USUBJID, VSTESTCD, VISIT, VSTPTNUM) %>%
#'   select(USUBJID, VSTESTCD, VISIT, VSTPTNUM, VSSTRESN, LASTFL)
#'
#' # Baseline (ABLFL) examples:
#'
#' input <- tribble(
#'   ~STUDYID,    ~USUBJID,  ~PARAMCD,     ~AVISIT,          ~ADT, ~AVAL,    ~DTYPE,
#'   "PILOT01",  "10-1187", "PARAM01",  "BASELINE",  "2021-04-27",  15.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM01",  "BASELINE",  "2021-04-25",  14.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM01",  "BASELINE",  "2021-04-23",  15.0, "AVERAGE",
#'   "PILOT01",  "10-1187", "PARAM01",    "WEEK 1",  "2021-04-27",  10.0, "AVERAGE",
#'   "PILOT01",  "10-1187", "PARAM01",    "WEEK 2",  "2021-04-30",  12.0,        NA,
#'   "PILOT01",  "10-1271", "PARAM01", "SCREENING",  "2021-04-27",  15.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM01",  "BASELINE",  "2021-04-25",  14.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM01",  "BASELINE",  "2021-04-23",  15.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM01",    "WEEK 1",  "2021-04-27",  10.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM01",    "WEEK 2",  "2021-04-30",  12.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM02", "SCREENING",  "2021-04-27",  15.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM02", "SCREENING",  "2021-04-25",  14.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM02", "SCREENING",  "2021-04-23",  15.0,        NA,
#'   "PILOT01",  "10-1271", "PARAM02",  "BASELINE",  "2021-04-27",  10.0, "AVERAGE",
#'   "PILOT01",  "10-1271", "PARAM02",    "WEEK 2",  "2021-04-30",  12.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM02", "SCREENING",  "2021-04-27",  15.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM02",  "BASELINE",  "2021-04-25",  14.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM02",    "WEEK 1",  "2021-04-23",  15.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM02",    "WEEK 1",  "2021-04-27",  10.0,        NA,
#'   "PILOT01",  "10-1187", "PARAM02",  "BASELINE",  "2021-04-30",  12.0,        NA
#' ) %>%
#'   mutate(
#'     ADT = as.Date(ADT)
#'   )
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
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AESEQ,        ~AEDECOD,     ~AEBODSYS,     ~AESEV, ~AESTDY,
#'   "PILOT01",    "AE", "10-1187",      1,    "IRRITATION",     "GENERAL",     "MILD",      44,
#'   "PILOT01",    "AE", "10-1271",      3,  "FIBRILLATION",     "CARDIAC", "MODERATE",      56,
#'   "PILOT01",    "AE", "10-1271",      4,       "CARDIAC",     "CARDIAC", "MODERATE",      57,
#'   "PILOT01",    "AE", "10-1271",      1,      "DYSPNOEA", "RESPIRATORY", "MODERATE",      56,
#'   "PILOT01",    "AE", "10-1271",      5, "HYPONATRAEMIA",  "METABOLISM", "MODERATE",      57,
#'   "PILOT01",    "AE", "10-1271",      2,    "MYOCARDIAL",     "CARDIAC",   "SEVERE",      56
#' )
#'
#' # Most severe AE first occurrence per patient
#' ae %>%
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
#' ae %>%
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
  assert_order_vars(order)
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
