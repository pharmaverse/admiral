#' Create a `event` Object
#'
#' The `event` object is used to define events as input for the
#' `derive_extreme_event()` and `derive_vars_extreme_event()` functions.
#'
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_extreme_event()`. If the argument is not specified, the input
#'   dataset (`dataset`) of `derive_extreme_event()` is used.
#'
#'   *Permitted Values*: a character scalar
#'
#' @param condition An unquoted condition for selecting the observations, which
#'   will contribute to the extreme event. If the condition contains summary
#'   functions like `all()`, they are evaluated for each by group separately.
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param mode If specified, the first or last observation with respect to `order` is
#'   selected for each by group.
#'
#'   *Permitted Values*: `"first"`, `"last"`, `NULL`
#'
#' @param order The specified variables or expressions are used to select the
#'   first or last observation if `mode` is specified.
#'
#'   `r roxygen_order_na_handling()`
#'
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @param set_values_to A named list returned by `exprs()` defining the variables
#'   to be set for the event, e.g. `exprs(PARAMCD = "WSP",
#'   PARAM  = "Worst Sleeping Problems")`. The values can be a symbol, a
#'   character string, a numeric value, `NA` or an expression.
#'
#'   *Permitted Values*: a named list of expressions, e.g., created by `exprs()`
#'
#' @param keep_source_vars Variables to keep from the source dataset
#'
#'   The specified variables are kept for the selected observations. The
#'   variables specified for `by_vars` (of `derive_extreme_event()`) and created
#'   by `set_values_to` are always kept.
#'
#'   *Permitted Values*: A list of expressions where each element is
#'   a symbol or a tidyselect expression, e.g., `exprs(VISIT, VISITNUM,
#'   starts_with("RS"))`.
#'
#' @param description Description of the event
#'
#'   The description does not affect the derivations where the event is used. It
#'   is intended for documentation only.
#'
#'   *Permitted Values*: a character scalar
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_extreme_event()], [derive_vars_extreme_event()], [event_joined()]
#'
#' @export
#'
#' @return An object of class `event`
event <- function(dataset_name = NULL,
                  condition = NULL,
                  mode = NULL,
                  order = NULL,
                  set_values_to = NULL,
                  keep_source_vars = NULL,
                  description = NULL) {
  out <- list(
    description = assert_character_scalar(description, optional = TRUE),
    dataset_name = assert_character_scalar(dataset_name, optional = TRUE),
    condition = assert_filter_cond(enexpr(condition), optional = TRUE),
    mode = assert_character_scalar(
      mode,
      values = c("first", "last"),
      case_sensitive = FALSE,
      optional = TRUE
    ),
    order = assert_expr_list(order, optional = TRUE),
    set_values_to = assert_expr_list(
      set_values_to,
      named = TRUE,
      optional = TRUE
    ),
    keep_source_vars = assert_expr_list(keep_source_vars, optional = TRUE)
  )
  class(out) <- c("event", "event_def", "source", "list")
  out
}

#' Create a `event_joined` Object
#'
#' @description
#'
#' The `event_joined` object is used to define events as input for the
#' `derive_extreme_event()` and `derive_vars_extreme_event()` functions.
#' This object should be used if the event does not depend on a single
#' observation of the source dataset but on multiple observations. For example,
#' if the event needs to be confirmed by a second observation of the source
#' dataset.
#'
#' The events are selected by calling `filter_joined()`. See its documentation
#' for more details.
#'
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_extreme_event()`. If the argument is not specified, the input
#'   dataset (`dataset`) of `derive_extreme_event()` is used.
#'
#'   *Permitted Values*: a character scalar
#'
#' @param condition An unquoted condition for selecting the observations, which
#'   will contribute to the extreme event.
#'
#'   The condition is applied to the joined dataset for selecting the confirmed
#'   observations. The condition can include summary functions like `all()` or
#'   `any()`. The joined dataset is grouped by the original observations. I.e.,
#'   the summary function are applied to all observations up to the confirmation
#'   observation. For example in the oncology setting when using this function
#'   for confirmed best overall response,  `condition = AVALC == "CR" &
#'   all(AVALC.join %in% c("CR", "NE")) & count_vals(var = AVALC.join, val =
#'   "NE") <= 1` selects observations with response "CR" and for all
#'   observations up to the confirmation observation the response is "CR" or
#'   "NE" and there is at most one "NE".
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param join_vars Variables to keep from joined dataset
#'
#'   The variables needed from the other observations should be specified for
#'   this parameter. The specified variables are added to the joined dataset
#'   with suffix ".join". For example to select all observations with `AVALC ==
#'   "Y"` and `AVALC == "Y"` for at least one subsequent visit `join_vars =
#'   exprs(AVALC, AVISITN)` and `condition = AVALC == "Y" & AVALC.join == "Y" &
#'   AVISITN < AVISITN.join` could be specified.
#'
#'   The `*.join` variables are not included in the output dataset.
#'
#'   *Permitted Values*: a named list of expressions, e.g., created by `exprs()`
#'
#' @param join_type Observations to keep after joining
#'
#'   The argument determines which of the joined observations are kept with
#'   respect to the original observation. For example, if `join_type =
#'   "after"` is specified all observations after the original observations are
#'   kept.
#'
#'   *Permitted Values:* `"before"`, `"after"`, `"all"`
#'
#' @param first_cond Condition for selecting range of data
#'
#'   `r lifecycle::badge("deprecated")`
#'
#'   This argument is *deprecated*, please use `first_cond_upper` instead.
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the subsequent observations, all
#'   observations are removed.
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param first_cond_lower Condition for selecting range of data (before)
#'
#'   If this argument is specified, the other observations are restricted from
#'   the first observation before the current observation where the specified
#'   condition is fulfilled up to the current observation. If the condition is
#'   not fulfilled for any of the other observations, no observations are
#'   considered, i.e., the observation is not flagged.
#'
#'   This parameter should be specified if `condition` contains summary
#'   functions which should not apply to all observations but only from a
#'   certain observation before the current observation up to the current
#'   observation.
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param first_cond_upper Condition for selecting range of data (after)
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the other observations, no
#'   observations are considered, i.e., the observation is not flagged.
#'
#'   This parameter should be specified if `condition` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment.
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param order If specified, the specified variables or expressions are used to
#'   select the first observation.
#'
#'   `r roxygen_order_na_handling()`
#'
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @inheritParams event
#'
#' @return An object of class `event_joined`
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_extreme_event()], [derive_vars_extreme_event()], [event()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#' library(lubridate)
#' # Derive confirmed best overall response (using event_joined())
#' # CR - complete response, PR - partial response, SD - stable disease
#' # NE - not evaluable, PD - progressive disease
#' adsl <- tribble(
#'   ~USUBJID, ~TRTSDTC,
#'   "1",      "2020-01-01",
#'   "2",      "2019-12-12",
#'   "3",      "2019-11-11",
#'   "4",      "2019-12-30",
#'   "5",      "2020-01-01",
#'   "6",      "2020-02-02",
#'   "7",      "2020-02-02",
#'   "8",      "2020-02-01"
#' ) %>%
#'   mutate(TRTSDT = ymd(TRTSDTC))
#'
#' adrs <- tribble(
#'   ~USUBJID, ~ADTC,        ~AVALC,
#'   "1",      "2020-01-01", "PR",
#'   "1",      "2020-02-01", "CR",
#'   "1",      "2020-02-16", "NE",
#'   "1",      "2020-03-01", "CR",
#'   "1",      "2020-04-01", "SD",
#'   "2",      "2020-01-01", "SD",
#'   "2",      "2020-02-01", "PR",
#'   "2",      "2020-03-01", "SD",
#'   "2",      "2020-03-13", "CR",
#'   "4",      "2020-01-01", "PR",
#'   "4",      "2020-03-01", "NE",
#'   "4",      "2020-04-01", "NE",
#'   "4",      "2020-05-01", "PR",
#'   "5",      "2020-01-01", "PR",
#'   "5",      "2020-01-10", "PR",
#'   "5",      "2020-01-20", "PR",
#'   "6",      "2020-02-06", "PR",
#'   "6",      "2020-02-16", "CR",
#'   "6",      "2020-03-30", "PR",
#'   "7",      "2020-02-06", "PR",
#'   "7",      "2020-02-16", "CR",
#'   "7",      "2020-04-01", "NE",
#'   "8",      "2020-02-16", "PD"
#' ) %>%
#'   mutate(
#'     ADT = ymd(ADTC),
#'     PARAMCD = "OVR",
#'     PARAM = "Overall Response by Investigator"
#'   ) %>%
#'   derive_vars_merged(
#'     dataset_add = adsl,
#'     by_vars = exprs(USUBJID),
#'     new_vars = exprs(TRTSDT)
#'   )
#'
#' derive_extreme_event(
#'   adrs,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(ADT),
#'   mode = "first",
#'   source_datasets = list(adsl = adsl),
#'   events = list(
#'     event_joined(
#'       description = paste(
#'         "CR needs to be confirmed by a second CR at least 28 days later",
#'         "at most one NE is acceptable between the two assessments"
#'       ),
#'       join_vars = exprs(AVALC, ADT),
#'       join_type = "after",
#'       first_cond_upper = AVALC.join == "CR" &
#'         ADT.join >= ADT + 28,
#'       condition = AVALC == "CR" &
#'         all(AVALC.join %in% c("CR", "NE")) &
#'         count_vals(var = AVALC.join, val = "NE") <= 1,
#'       set_values_to = exprs(
#'         AVALC = "CR"
#'       )
#'     ),
#'     event_joined(
#'       description = paste(
#'         "PR needs to be confirmed by a second CR or PR at least 28 days later,",
#'         "at most one NE is acceptable between the two assessments"
#'       ),
#'       join_vars = exprs(AVALC, ADT),
#'       join_type = "after",
#'       first_cond_upper = AVALC.join %in% c("CR", "PR") &
#'         ADT.join >= ADT + 28,
#'       condition = AVALC == "PR" &
#'         all(AVALC.join %in% c("CR", "PR", "NE")) &
#'         count_vals(var = AVALC.join, val = "NE") <= 1,
#'       set_values_to = exprs(
#'         AVALC = "PR"
#'       )
#'     ),
#'     event(
#'       description = paste(
#'         "CR, PR, or SD are considered as SD if occurring at least 28",
#'         "after treatment start"
#'       ),
#'       condition = AVALC %in% c("CR", "PR", "SD") & ADT >= TRTSDT + 28,
#'       set_values_to = exprs(
#'         AVALC = "SD"
#'       )
#'     ),
#'     event(
#'       condition = AVALC == "PD",
#'       set_values_to = exprs(
#'         AVALC = "PD"
#'       )
#'     ),
#'     event(
#'       condition = AVALC %in% c("CR", "PR", "SD", "NE"),
#'       set_values_to = exprs(
#'         AVALC = "NE"
#'       )
#'     ),
#'     event(
#'       description = "set response to MISSING for patients without records in ADRS",
#'       dataset_name = "adsl",
#'       condition = TRUE,
#'       set_values_to = exprs(
#'         AVALC = "MISSING"
#'       ),
#'       keep_source_vars = exprs(TRTSDT)
#'     )
#'   ),
#'   set_values_to = exprs(
#'     PARAMCD = "CBOR",
#'     PARAM = "Best Confirmed Overall Response by Investigator"
#'   )
#' ) %>%
#'   filter(PARAMCD == "CBOR")
event_joined <- function(dataset_name = NULL,
                         condition,
                         order = NULL,
                         join_vars,
                         join_type,
                         first_cond = NULL,
                         first_cond_lower = NULL,
                         first_cond_upper = NULL,
                         set_values_to = NULL,
                         keep_source_vars = NULL,
                         description = NULL) {
  if (!missing(first_cond)) {
    deprecate_stop(
      "1.1.0",
      "event_joined(first_cond=)",
      "event_joined(first_cond_upper=)"
    )
    first_cond_upper <- assert_filter_cond(enexpr(first_cond), optional = TRUE)
  } else {
    first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  }

  out <- list(
    description = assert_character_scalar(description, optional = TRUE),
    dataset_name = assert_character_scalar(dataset_name, optional = TRUE),
    condition = assert_filter_cond(enexpr(condition), optional = TRUE),
    order = assert_expr_list(order, optional = TRUE),
    join_vars = assert_vars(join_vars),
    join_type = assert_character_scalar(
      join_type,
      values = c("before", "after", "all"),
      case_sensitive = FALSE
    ),
    first_cond_lower = assert_filter_cond(enexpr(first_cond_lower), optional = TRUE),
    first_cond_upper = first_cond_upper,
    set_values_to = assert_expr_list(
      set_values_to,
      named = TRUE,
      optional = TRUE
    ),
    keep_source_vars = assert_expr_list(keep_source_vars, optional = TRUE)
  )
  class(out) <- c("event_joined", "event_def", "source", "list")
  out
}
