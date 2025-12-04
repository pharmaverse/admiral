#' Add a Variable Flagging the First or Last Observation Within Each By Group
#'
#' Add a variable flagging the first or last observation within each by group
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param order Sort order
#'
#'   The first or last observation is determined with respect to the specified
#'   order.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [var_list]
#'
#' @param new_var Variable to add
#'
#'   The specified variable is added to the output dataset. It is set to the value
#'   set in `true_value` for the first or last observation (depending on the mode) of each by group.
#'
#' @permitted [var]
#'
#' @param mode Flag mode
#'
#'   Determines of the first or last observation is flagged.
#'
#' @permitted [mode]
#'
#' @param true_value True value
#'
#'   The value for the specified variable `new_var`, applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#' @permitted [char_scalar]
#'
#' @param false_value False value
#'
#'   The value for the specified variable `new_var`, NOT applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#' @permitted [char_scalar]
#'
#' @param flag_all Flag setting
#'
#'   A logical value where if set to `TRUE`, all records are flagged
#'   and no error or warning is issued if the first or last record is not unique.
#'
#' @permitted [boolean]
#'
#' @param by_vars Grouping variables
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#' @permitted [msg_type]
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter), `new_var` is set to `"Y"` for the first or last observation
#'   (with respect to the order specified for the `order` parameter and the flag mode
#'   specified for the `mode` parameter). In the case where the user wants to flag multiple records
#'   of a grouping, for example records that all happen on the same visit and time, the argument
#'   `flag_all` can be set to `TRUE`.
#'   Otherwise, `new_var` is set to `NA`. Thus, the direction of "worst" is considered fixed for
#'   all parameters in the dataset depending on the `order` and the `mode`, i.e. for every
#'   parameter the first or last record will be flagged across the whole dataset.
#'
#'
#' @return The input dataset with the new flag variable added
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examplesx
#'
#' @caption Data setup
#'
#' @info The following examples use the `ADVS` and `ADAE` datasets below as a basis.
#'
#' @code
#' library(tibble, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD,    ~AVISIT,          ~ADT, ~AVAL,
#'   "1015",   "TEMP",   "BASELINE",  "2021-04-27",  38.0,
#'   "1015",   "TEMP",   "BASELINE",  "2021-04-25",  39.0,
#'   "1015",   "TEMP",   "WEEK 2",    "2021-05-10",  37.5,
#'   "1015",   "WEIGHT", "SCREENING", "2021-04-19",  81.2,
#'   "1015",   "WEIGHT", "BASELINE",  "2021-04-25",  82.7,
#'   "1015",   "WEIGHT", "BASELINE",  "2021-04-27",  84.0,
#'   "1015",   "WEIGHT", "WEEK 2",    "2021-05-09",  82.5,
#'   "1023",   "TEMP",   "SCREENING", "2021-04-27",  38.0,
#'   "1023",   "TEMP",   "BASELINE",  "2021-04-28",  37.5,
#'   "1023",   "TEMP",   "BASELINE",  "2021-04-29",  37.5,
#'   "1023",   "TEMP",   "WEEK 1",    "2021-05-03",  37.0,
#'   "1023",   "WEIGHT", "SCREENING", "2021-04-27",  69.6,
#'   "1023",   "WEIGHT", "BASELINE",  "2021-04-29",  67.2,
#'   "1023",   "WEIGHT", "WEEK 1",    "2021-05-02",  65.9
#' ) %>%
#' mutate(
#'   STUDYID = "AB123",
#'   ADT = ymd(ADT)
#' )
#'
#' adae <- tribble(
#'   ~USUBJID,         ~AEBODSYS,    ~AEDECOD,   ~AESEV, ~AESTDY, ~AESEQ,
#'   "1015", "GENERAL DISORDERS",  "ERYTHEMA",   "MILD",       2,      1,
#'   "1015", "GENERAL DISORDERS",  "PRURITUS",   "MILD",       2,      2,
#'   "1015",      "GI DISORDERS", "DIARRHOEA",   "MILD",       8,      3,
#'   "1023", "CARDIAC DISORDERS",  "AV BLOCK",   "MILD",      22,      4,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       3,      1,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA", "SEVERE",       5,      2,
#'   "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       8,      3
#' ) %>%
#' mutate(STUDYID = "AB123")
#'
#' @caption Flagging the first/last observation within a by group (`order`, `mode`)
#'
#' @info A new variable is added for each subject to flag the last observation
#'   within a by group. Within each by group (specified by `by_vars`), the
#'   `order = exprs(ADT)` argument specifies we wish to sort the records by analysis
#'   date and then select the last one (`mode = "last"`). The name of the new
#'   variable is passed through the `new_var = LASTFL` call.
#'
#' @code
#' advs %>%
#'   derive_var_extreme_flag(
#'     by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'     order = exprs(ADT),
#'     new_var = LASTFL,
#'     mode = "last",
#'   ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())
#'
#' @info Note here that a similar `FIRSTFL` variable could instead be derived
#'    simply by switching to `mode = "first"`. Alternatively, we could make use
#'    of `desc()` within the sorting specified by `order`:
#'
#' @code
#' advs %>%
#'   derive_var_extreme_flag(
#'     by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'     order = exprs(desc(ADT)),
#'     new_var = FIRSTFL,
#'     mode = "last",
#'   ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())
#'
#' @caption Modifying the flag values (`true_value`, `false_value`)
#'
#' @info The previous example is now enhanced with custom values for the flag
#'   entries. Records which are flagged are filled with the contents of
#'   `true_value` and those which are not are filled with the contents of
#'   `false_value`. Note that these are normally preset to `"Y"` and `NA`,
#'   which is why they were not specified in the example above.
#'
#' @code
#' advs %>%
#'   derive_var_extreme_flag(
#'     by_vars = exprs(STUDYID, USUBJID, PARAMCD),
#'     order = exprs(ADT),
#'     new_var = LASTFL,
#'     mode = "last",
#'     true_value = "Yes",
#'     false_value = "No",
#'   ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())
#'
#' @caption Creating temporary variables for sorting (`check_type`)
#'
#' @info In this example we wish to flag the first occurrence of the most
#'   severe AE within each subject. To ensure correct sorting of the
#'   severity values, `AESEV` must be pre-processed into a numeric variable
#'   `TEMP_AESEVN` which can then be passed inside `order`. Once again,
#'   to ensure we only flag the *first* occurrence, we specify `AESTDY` and
#'   `AESEQ` inside `order` as well.
#'
#' @code
#' adae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCIFL,
#'     by_vars = exprs(STUDYID, USUBJID),
#'     order = exprs(TEMP_AESEVN, AESTDY, AESEQ),
#'     mode = "first",
#'     check_type = "warning"
#'   ) %>%
#'   arrange(STUDYID, USUBJID, AESTDY, AESEQ) %>%
#'   select(STUDYID, USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
#'
#' @info Note here that the presence of `AESEQ` as a sorting variable inside
#'   the `order` argument ensures that the combination of `by_vars` and
#'   `order` indexes unique records in the dataset. If this had been omitted,
#'   the choice of `check_type = "warning"` would have ensured that
#'   `derive_var_extreme_flag()` would throw a warning due to perceived
#'   duplicate records (in this case, the first two AEs for subject `"1015"`).
#'   If no sorting variables exist, or if these duplicates are acceptable, then
#'   the user can silence the warning with `check_type = "none"`. Alternatively,
#'   the warning can be upgraded to an error with `check_type = "error"`.
#'
#' @caption Flagging all records if multiple are identified (`flag_all`)
#'
#' @info Revisiting the above example, if we instead wish to flag *all* AEs
#'   of the highest severity occurring on the earliest date, then we can use
#'   `flag_all = TRUE`. Note that we now also omit `AESEQ` from the `order`
#'   argument because we do not need to differentiate between two AEs occurring
#'   on the same day (e.g. for subject `"1015"`) as they are both flagged.
#'
#' @code
#' adae %>%
#'   mutate(
#'     TEMP_AESEVN =
#'       as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
#'   ) %>%
#'   derive_var_extreme_flag(
#'     new_var = AOCCIFL,
#'     by_vars = exprs(STUDYID, USUBJID),
#'     order = exprs(TEMP_AESEVN, AESTDY),
#'     mode = "first",
#'     flag_all = TRUE
#'   ) %>%
#'   arrange(STUDYID, USUBJID, AESTDY, AESEQ) %>%
#'   select(STUDYID, USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
#'
#' @caption Deriving a baseline flag
#'
#' @info `derive_var_extreme_flag()` is very often used to derive the baseline
#'   flag `ABLFL`, so the following section contains various examples of this
#'   in action for the `ADVS` dataset. Note that for these derivations it is
#'   often convenient to leverage higher order functions such as
#'   `restrict_derivation()` and `slice_derivation()`. Please read the
#'   [Higher Order Functions](https://pharmaverse.github.io/admiral/articles/higher_order.html)
#'   vignette, as well as their specific reference pages, to learn more.
#'
#'   To set the baseline flag for the last observation among those where
#'   `AVISIT = "BASELINE"`, we can use a similar call to the examples above but
#'   wrapping inside of `restrict_derivation()` and making use of the `filter`
#'   argument.
#'
#' @code
#' restrict_derivation(
#'   advs,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(ADT),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE"
#' ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())
#'
#' @info Alternatively, to set baseline as the lowest observation among
#'   those where `AVISIT = "BASELINE"` (selecting the latest if there are
#'   multiple) we can modify the `order` argument, ensuring to sort by
#'   descending `AVAL`  before `ADT`. Note here the synergy between `desc()`
#'   and `mode`, because `mode = "last"` applies to both the ordering
#'   variables `AVAL` and `ADT` and so we need to reverse only the ordering
#'   of the former to ensure that the lowest value is selected but also that
#'   the latest one among multiple is preferred. This is relevant for
#'   subject `"1023"`'s temperature records.
#'
#' @code
#' restrict_derivation(
#'   advs,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     order = exprs(desc(AVAL), ADT),
#'     new_var = ABLFL,
#'     mode = "last"
#'   ),
#'   filter = AVISIT == "BASELINE"
#' ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())
#'
#' @info In practice, baseline-setting may vary on a parameter by
#'   parameter basis, in which case `slice_derivation()` could be
#'   used in place of `restrict_derivation()`. In the example below, we
#'   set the baseline flag as follows: for temperature records, as the
#'   lowest value recorded at a baseline visit; for weight records,
#'   as the highest value recorded at a baseline visit. In both cases,
#'   we again select the latest observation if there are multiple.
#'
#' @code
#' slice_derivation(
#'   advs,
#'   derivation = derive_var_extreme_flag,
#'   args = params(
#'     by_vars = exprs(USUBJID, PARAMCD),
#'     mode = "last",
#'     new_var = ABLFL,
#'   ),
#'   derivation_slice(
#'     filter = AVISIT == "BASELINE" & PARAMCD == "TEMP",
#'     args = params(order = exprs(desc(AVAL), ADT))
#'   ),
#'   derivation_slice(
#'     filter = AVISIT == "BASELINE" & PARAMCD == "WEIGHT",
#'     args = params(order = exprs(AVAL, ADT))
#'   )
#' ) %>%
#'   arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
#'   select(STUDYID, everything())

derive_var_extreme_flag <- function(dataset,
                                    by_vars,
                                    order,
                                    new_var,
                                    mode,
                                    true_value = "Y",
                                    false_value = NA_character_,
                                    flag_all = FALSE,
                                    check_type = "warning") {
  new_var <- assert_symbol(enexpr(new_var))
  assert_vars(by_vars)
  assert_expr_list(order)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, !!!extract_vars(order)))
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  assert_atomic_vector(true_value, optional = TRUE)
  assert_atomic_vector(false_value, optional = TRUE)
  flag_all <- assert_logical_scalar(flag_all)
  check_type <- assert_character_scalar(
    check_type,
    values = c("none", "warning", "error"),
    case_sensitive = FALSE
  )

  # Create flag
  if (flag_all) {
    check_type <- "none"
  }

  # Create observation number to identify the extreme record
  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr,
      order = order,
      by_vars = by_vars,
      check_type = check_type
    )

  if (mode == "first") {
    data <- data %>%
      mutate(!!new_var := if_else(!!tmp_obs_nr == 1, true_value, false_value))
  } else {
    data <- data %>%
      group_by(!!!by_vars) %>%
      mutate(!!new_var := if_else(!!tmp_obs_nr == n(), true_value, false_value)) %>%
      ungroup()
  }

  if (flag_all) {
    flag_direction <- if_else(mode == "first", "down", "up")
    data <- data %>%
      group_by(!!!by_vars, !!!order) %>%
      fill(!!new_var, .direction = flag_direction) %>%
      ungroup()
  }


  # Remove temporary variable
  data %>%
    remove_tmp_vars()
}
