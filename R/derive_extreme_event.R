#' Add the Worst or Best Observation for Each By Group as New Records
#'
#' Add the first available record from `events` for each by group as new
#' records, all variables of the selected observation are kept. It can be used
#' for selecting the extreme observation from a series of user-defined events.
#' This distinguish `derive_extreme_event()` from `derive_extreme_records()`,
#' where extreme records are derived based on certain order of existing
#' variables.
#'
#' @param order Sort order
#'
#'   If a particular event from `events` has more than one observation, within
#'   the event and by group, the records are ordered by the specified order.
#'
#'   *Permitted Values:* list of variables or `desc(<variable>)` function calls
#'   created by `exprs()`, e.g., `exprs(ADT, desc(AVAL))`
#'
#'
#' @param mode Selection mode (first or last)
#'
#'   If a particular event from `events` has more than one observation,
#'   "first"/"last" is to select the first/last record of this type of events
#'   sorting by `order`.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param events Conditions and new values defining events
#'
#'   A list of `event()` objects is expected. Only observations listed in the
#'   `events` are considered for deriving extreme event. If multiple records
#'   meet the filter `condition`, take the first record sorted by `order`.
#'
#' @inheritParams filter_extreme
#' @inheritParams derive_summary_records
#'
#' @details
#'   1. Construct a dataset based on `events`: apply the filter `condition` and
#'   `set_values_to` to the input dataset.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is selected.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the selected observations.
#'   1. The observations are added to input dataset.
#'
#'
#' @return The input dataset with the best or worst observation of each by group
#'   added as new observations.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adqs <- tribble(
#'   ~USUBJID, ~PARAMCD,       ~AVALC,        ~ADY,
#'   "1",      "NO SLEEP",     "N",              1,
#'   "1",      "WAKE UP",      "N",              2,
#'   "1",      "FALL ASLEEP",  "N",              3,
#'   "2",      "NO SLEEP",     "N",              1,
#'   "2",      "WAKE UP",      "Y",              2,
#'   "2",      "WAKE UP",      "Y",              3,
#'   "2",      "FALL ASLEEP",  "N",              4,
#'   "3",      "NO SLEEP",     NA_character_,    1
#' )
#'
#' # Add a new record for each USUBJID storing the the worst sleeping problem.
#' derive_extreme_event(
#'   adqs,
#'   by_vars = exprs(USUBJID),
#'   events = list(
#'     event(
#'       condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
#'     ),
#'     event(
#'       condition = PARAMCD == "WAKE UP" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "Waking up more than three times", AVAL = 2)
#'     ),
#'     event(
#'       condition = PARAMCD == "FALL ASLEEP" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "More than 30 mins to fall asleep", AVAL = 3)
#'     ),
#'     event(
#'       condition = all(AVALC == "N"),
#'       set_values_to = exprs(
#'         AVALC = "No sleeping problems", AVAL = 4
#'       )
#'     ),
#'     event(
#'       condition = TRUE,
#'       set_values_to = exprs(AVALC = "Missing", AVAL = 99)
#'     )
#'   ),
#'   order = exprs(ADY),
#'   mode = "last",
#'   set_values_to = exprs(
#'     PARAMCD = "WSP",
#'     PARAM = "Worst Sleeping Problems"
#'   )
#' )
derive_extreme_event <- function(dataset,
                                 by_vars = NULL,
                                 events,
                                 order,
                                 mode,
                                 check_type = "warning",
                                 set_values_to) {
  # Check input parameters
  assert_vars(by_vars, optional = TRUE)
  assert_list_of(events, "event")
  assert_order_vars(order)
  assert_data_frame(
    dataset,
    required_vars = by_vars
  )
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_varval_list(set_values_to)

  # Create new observations
  ## Create a dataset (selected_records) from `events`
  condition_ls <- map(events, "condition")
  set_values_to_ls <- map(events, "set_values_to")
  event_order <- map(seq_len(length(events)), function(x) x)
  tmp_event_no <- get_new_tmp_var(dataset, prefix = "tmp_event_no")

  selected_records_ls <- purrr::pmap(
    list(condition_ls, set_values_to_ls, event_order),
    function(x, y, z) {
      dataset %>%
        group_by(!!!by_vars) %>%
        filter(!!enexpr(x)) %>%
        mutate(!!!y, !!tmp_event_no := z) %>%
        ungroup()
    }
  )
  selected_records <- bind_rows(selected_records_ls)

  ## tmp obs number within by_vars and a type of event
  tmp_obs <- get_new_tmp_var(selected_records)
  selected_records_obs <- selected_records %>%
    derive_var_obs_number(
      new_var = !!tmp_obs,
      order = order,
      by_vars = expr_c(by_vars, expr(!!tmp_event_no)),
      check_type = check_type
    )

  ## filter_extreme
  if (mode == "first") {
    tmp_obs_expr <- expr(!!tmp_obs)
  } else {
    tmp_obs_expr <- expr(desc(!!tmp_obs))
  }
  new_obs <- selected_records_obs %>%
    filter_extreme(
      by_vars = by_vars,
      order = expr_c(expr(!!tmp_event_no), tmp_obs_expr),
      mode = "first",
      check_type = check_type
    ) %>%
    mutate(!!!set_values_to) %>%
    select(-!!tmp_event_no, -!!tmp_obs)

  # Create output dataset
  bind_rows(dataset, new_obs)
}

#' Create a `event` Object
#'
#' The `event` object is used to define events as input for the
#' `derive_extreme_event()` function.
#'
#' @param condition An unquoted condition for selecting the observations, which
#'   will contribute to the extreme event.
#'
#' @param set_values_to A named list returned by `exprs()` defining the variables
#'   to be set for the extreme answer, e.g. `exprs(PARAMCD = "WSP",
#'   PARAM  = "Worst Sleeping Problems"`. The values must be a symbol, a
#'   character string, a numeric value, or `NA`.
#'
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_extreme_event()]
#'
#' @export
#'
#' @return An object of class `event`
event <- function(condition,
                  set_values_to = NULL) {
  out <- list(
    condition = assert_filter_cond(enexpr(condition), optional = TRUE),
    set_values_to = assert_varval_list(
      set_values_to,
      optional = TRUE
    )
  )
  class(out) <- c("event", "source", "list")
  out
}
