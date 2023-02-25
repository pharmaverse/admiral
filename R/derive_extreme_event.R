#' Add the Worst or Best Observation for Each By Group as New Records
#'
#' Add the worst or best observation for each by group as new records, all
#' variables of the selected observation are kept. It can be used for selecting
#' the extreme observation from a series of user-defined events. This
#' distinguish `derive_extreme_event()` from `derive_extreme_records()`, where
#' extreme records are derived based on certain order of existing variables.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first observation of each by group is added
#'   to the input dataset. If `"last"` is specified, the last observation of
#'   each by group is added to the input dataset.
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
#'   ~USUBJID, ~PARAMCD,       ~AVALC,
#'   "1",      "NO SLEEP",     "N",
#'   "1",      "WAKE UP",      "N",
#'   "1",      "FALL ASLEEP",  "N",
#'   "2",      "NO SLEEP",     "N",
#'   "2",      "WAKE UP",      "Y",
#'   "2",      "FALL ASLEEP",  "N",
#'   "3",      "NO SLEEP",     NA_character_
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
#'   order = exprs(AVAL),
#'   mode = "first",
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
  condition_ls <- map(events, 1)
  set_values_to_ls <- map(events, 2)
  selected_records_ls <- map2(
    condition_ls, set_values_to_ls,
    function(x, y) {
      dataset %>%
        group_by(!!!by_vars) %>%
        filter(!!enexpr(x)) %>%
        mutate(!!!y) %>%
        arrange(!!!order, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup()
    }
  )
  selected_records <- bind_rows(selected_records_ls)
  ## filter_extreme for selected_records
  new_obs <- selected_records %>%
    filter_extreme(
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    ) %>%
    mutate(!!!set_values_to)

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
