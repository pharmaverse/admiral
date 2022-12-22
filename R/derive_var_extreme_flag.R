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
#' @param filter Deprecated, please use `restrict_derivation()` instead (see examples).
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
#' @author Stefan Bundfuss
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
#' library(admiral.test)
#' data("admiral_vs")
#'
#' # Flag last value for each patient, test, and visit, baseline observations are ignored
#' admiral_vs %>%
#'   restrict_derivation(
#'     derivation = derive_var_extreme_flag,
#'     args = params(
#'       by_vars = vars(USUBJID, VSTESTCD, VISIT),
#'       order = vars(VSTPTNUM),
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
#'     by_vars = vars(USUBJID, PARAMCD),
#'     order = vars(ADT),
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
#'     by_vars = vars(USUBJID, PARAMCD),
#'     order = vars(AVAL, ADT),
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
#'     by_vars = vars(USUBJID, PARAMCD),
#'     order = vars(desc(AVAL), ADT),
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
#'     by_vars = vars(USUBJID, PARAMCD),
#'     order = vars(ADT, desc(AVAL)),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE" & DTYPE == "AVERAGE"
#' )
#'
#' # OCCURDS Examples
#' data("admiral_ae")
#'
#' # Most severe AE first occurrence per patient
#' admiral_ae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCIFL,
#'     by_vars = vars(USUBJID),
#'     order = vars(TEMP_AESEVN, AESTDY, AESEQ),
#'     mode = "first"
#'   ) %>%
#'   arrange(USUBJID, AESTDY, AESEQ) %>%
#'   select(USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
#'
#' # Most severe AE first occurrence per patient per body system
#' admiral_ae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCSIFL,
#'     by_vars = vars(USUBJID, AEBODSYS),
#'     order = vars(TEMP_AESEVN, AESTDY, AESEQ),
#'     mode = "first"
#'   ) %>%
#'   arrange(USUBJID, AESTDY, AESEQ) %>%
#'   select(USUBJID, AEBODSYS, AESEV, AESTDY, AOCCSIFL)
derive_var_extreme_flag <- function(dataset,
                                    by_vars,
                                    order,
                                    new_var,
                                    mode,
                                    filter = deprecated(),
                                    check_type = "warning") {
  if (!missing(filter)) {
    deprecate_stop("0.7.0", "derive_var_extreme_flag(filter = )", "restrict_derivation(filter = )")
  }

  new_var <- assert_symbol(enquo(new_var))
  assert_vars(by_vars)
  assert_order_vars(order)
  assert_data_frame(dataset, required_vars = vars(!!!by_vars, !!!extract_vars(order)))
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
#' is labelled in the `new_var` column as `"Y"`,
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
#' @author Ondrej Slama
#'
#' @return The input dataset with the new flag variable added.
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' input <- tribble(
#'   ~STUDYID, ~USUBJID,  ~PARAMCD,     ~AVISIT,                  ~ADT, ~AVAL,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-27"),  15.0,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-25"),  14.0,
#'   "TEST01",  "PAT01", "PARAM01",  "BASELINE", as.Date("2021-04-23"),  15.0,
#'   "TEST01",  "PAT01", "PARAM01",    "WEEK 1", as.Date("2021-04-27"),  10.0,
#'   "TEST01",  "PAT01", "PARAM01",    "WEEK 2", as.Date("2021-04-30"),  12.0,
#'   "TEST01",  "PAT02", "PARAM01", "SCREENING", as.Date("2021-04-27"),  15.0,
#'   "TEST01",  "PAT02", "PARAM01",  "BASELINE", as.Date("2021-04-25"),  14.0,
#'   "TEST01",  "PAT02", "PARAM01",  "BASELINE", as.Date("2021-04-23"),  15.0,
#'   "TEST01",  "PAT02", "PARAM01",    "WEEK 1", as.Date("2021-04-27"),  10.0,
#'   "TEST01",  "PAT02", "PARAM01",    "WEEK 2", as.Date("2021-04-30"),  12.0,
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-27"),  15.0,
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-25"),  14.0,
#'   "TEST01",  "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-23"),  15.0,
#'   "TEST01",  "PAT01", "PARAM02",  "BASELINE", as.Date("2021-04-27"),  10.0,
#'   "TEST01",  "PAT01", "PARAM02",    "WEEK 2", as.Date("2021-04-30"),  12.0,
#'   "TEST01",  "PAT02", "PARAM02", "SCREENING", as.Date("2021-04-27"),  15.0,
#'   "TEST01",  "PAT02", "PARAM02",  "BASELINE", as.Date("2021-04-25"),  14.0,
#'   "TEST01",  "PAT02", "PARAM02",    "WEEK 1", as.Date("2021-04-23"),  15.0,
#'   "TEST01",  "PAT02", "PARAM02",    "WEEK 1", as.Date("2021-04-27"),  10.0,
#'   "TEST01",  "PAT02", "PARAM02",  "BASELINE", as.Date("2021-04-30"),  12.0,
#'   "TEST01",  "PAT02", "PARAM03", "SCREENING", as.Date("2021-04-27"),  15.0,
#'   "TEST01",  "PAT02", "PARAM03",  "BASELINE", as.Date("2021-04-25"),  14.0,
#'   "TEST01",  "PAT02", "PARAM03",    "WEEK 1", as.Date("2021-04-23"),  15.0,
#'   "TEST01",  "PAT02", "PARAM03",    "WEEK 1", as.Date("2021-04-27"),  10.0,
#'   "TEST01",  "PAT02", "PARAM03",  "BASELINE", as.Date("2021-04-30"),  12.0
#' )
#'
#' derive_var_worst_flag(
#'   input,
#'   by_vars = vars(USUBJID, PARAMCD, AVISIT),
#'   order = vars(desc(ADT)),
#'   new_var = WORSTFL,
#'   param_var = PARAMCD,
#'   analysis_var = AVAL,
#'   worst_high = c("PARAM01", "PARAM03"),
#'   worst_low = "PARAM02"
#' )
#' \dontrun{
#' # example with ADVS
#' restrict_derivation(
#'   advs,
#'   derivation = derive_var_worst_flag,
#'   args = params(
#'     by_vars = vars(USUBJID, PARAMCD, AVISIT),
#'     order = vars(ADT, ATPTN),
#'     new_var = WORSTFL,
#'     param_var = PARAMCD,
#'     analysis_var = AVAL,
#'     worst_high = c("SYSBP", "DIABP"),
#'     worst_low = "RESP"
#'   ),
#'   filter = !is.na(AVISIT) & !is.na(AVAL)
#' )
#' }
#'
derive_var_worst_flag <- function(dataset,
                                  by_vars,
                                  order,
                                  new_var,
                                  param_var,
                                  analysis_var,
                                  worst_high,
                                  worst_low,
                                  filter = deprecated(),
                                  check_type = "warning") {
  if (!missing(filter)) {
    deprecate_stop("0.7.0", "derive_var_worst_flag(filter = )", "restrict_derivation(filter = )")
  }

  # perform argument checks
  new_var <- assert_symbol(enquo(new_var))
  param_var <- assert_symbol(enquo(param_var))
  analysis_var <- assert_symbol(enquo(analysis_var))
  assert_vars(by_vars)
  assert_order_vars(order)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, extract_vars(order), param_var, analysis_var)
  )
  assert_character_vector(worst_high)
  assert_character_vector(worst_low)

  # additional checks for worstflag - parameters overlap
  if (length(intersect(worst_high, worst_low)) > 0) {
    err_msg <- paste(
      "The following parameter(-s) are both assigned to `worst_high` and `worst_low` flags:",
      paste0(intersect(worst_high, worst_low), collapse = ", ")
    )
    abort(err_msg)
  }

  # additional checks for worstflag - parameters not available
  param_var_str <- as_string(quo_get_expr(param_var))
  if (length(worst_high) > 0 &&
    !all(worst_high %in% dataset[[param_var_str]])) {
    err_msg <- paste0(
      "The following parameter(-s) in `worst_high` are not available in column ",
      param_var_str,
      ": ",
      paste0(worst_high[!worst_high %in% dataset[[param_var_str]]], collapse = ", ")
    )
    abort(err_msg)
  }

  # additional checks for worstflag - parameters not available
  if (length(worst_low) > 0 &&
    !all(worst_low %in% dataset[[param_var_str]])) {
    err_msg <- paste0(
      "The following parameter(-s) in `worst_low` are not available in column ",
      param_var_str,
      ": ",
      paste0(worst_low[!worst_low %in% dataset[[param_var_str]]], collapse = ", ")
    )
    abort(err_msg)
  }

  # derive worst-flag
  bind_rows(
    derive_var_extreme_flag(
      dataset = dplyr::filter(dataset, !!param_var %in% worst_low),
      by_vars = by_vars,
      order = quo_c(analysis_var, order),
      new_var = !!new_var,
      mode = "first",
      check_type = check_type
    ),
    derive_var_extreme_flag(
      dataset = dplyr::filter(dataset, !!param_var %in% worst_high),
      by_vars = by_vars,
      order = quo_c(quo(desc(!!quo_get_expr(analysis_var))), order),
      new_var = !!new_var,
      mode = "first",
      check_type = check_type
    ),
    dplyr::filter(dataset, !(!!param_var %in% c(worst_low, worst_high)))
  )
}
