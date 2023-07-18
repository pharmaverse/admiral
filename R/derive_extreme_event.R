#' Add the Worst or Best Observation for Each By Group as New Records
#'
#' Add the first available record from `events` for each by group as new
#' records, all variables of the selected observation are kept. It can be used
#' for selecting the extreme observation from a series of user-defined events.
#' This distinguish `derive_extreme_event()` from `derive_extreme_records()`,
#' where extreme records are derived based on certain order of existing
#' variables.
#'
#' @param events Conditions and new values defining events
#'
#'   A list of `event()` or `event_joined()` objects is expected. Only
#'   observations listed in the `events` are considered for deriving extreme
#'   event. If multiple records meet the filter `condition`, take the first
#'   record sorted by `order`. The data is grouped by `by_vars`, i.e., summary
#'   functions like `all()` or `any()` can be used in `condition`.
#'
#'   For `event_joined()` events the observations are selected by calling
#'   `filter_joined`. The `condition` field is passed to the `filter` argument.
#'
#' @param order Sort order
#'
#'   If a particular event from `events` has more than one observation, within
#'   the event and by group, the records are ordered by the specified order.
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If a particular event from `events` has more than one observation,
#'   "first"/"last" is to select the first/last record of this type of events
#'   sorting by `order`.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of `event()`
#'   and `event_joined()` refers to the dataset provided in the list.
#'
#' @param keep_vars_source Variables to keep from the source dataset
#'
#'   For each event the specified variables are kept from the selected
#'   observations. The variables specified for `by_vars` and created by
#'   `set_values_to` are always kept.
#'
#'   *Permitted Values*: A list of expressions where each element is
#'   a symbol or a tidyselect expression, e.g., `exprs(VISIT, VISITNUM,
#'   starts_with("RS"))`.
#'
#' @inheritParams filter_extreme
#' @inheritParams derive_summary_records
#'
#' @details
#'   1. For each event select the observations to consider:
#'
#'       1. If the event is of class `event`, the observations of the source dataset
#'       are restricted by `condition` and then the first or last (`mode`)
#'       observation per by group (`by_vars`) is selected.
#'
#'           If the event is of class `event_joined`, `filter_joined()` is called to
#'           select the observations.
#'
#'       1. The variables specified by the `set_values_to` field of the event
#'       are added to the selected observations.
#'       1. Only the variables specified for the `keep_vars_source` field of the
#'       event, and the by variables (`by_vars`) and the variables created by
#'       `set_values_to` are kept.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first event is selected. If there is more than one
#'   observation per event the first or last observation (with respect to the
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
                                 source_datasets = NULL,
                                 check_type = "warning",
                                 set_values_to,
                                 keep_vars_source = exprs(everything())) {
  # Check input parameters
  assert_vars(by_vars, optional = TRUE)
  assert_list_of(events, "event_def")
  assert_expr_list(order)
  assert_data_frame(
    dataset,
    required_vars = by_vars
  )
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  assert_list_of(source_datasets, "data.frame")
  source_names <- names(source_datasets)
  events_to_check <- events[map_lgl(events, ~ !is.null(.x$dataset_name))]
  if (length(events_to_check) > 0) {
    assert_list_element(
      list = events_to_check,
      element = "dataset_name",
      condition = dataset_name %in% source_names,
      source_names = source_names,
      message_text = paste0(
        "The dataset names must be included in the list specified for the ",
        "`source_datasets` parameter.\n",
        "Following names were provided by `source_datasets`:\n",
        enumerate(source_names, quote_fun = squote)
      )
    )
  }

  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_varval_list(set_values_to)
  keep_vars_source <- assert_expr_list(keep_vars_source)

  # Create new observations
  ## Create a dataset (selected_records) from `events`
  event_index <- as.list(seq_along(events))
  tmp_event_no <- get_new_tmp_var(dataset, prefix = "tmp_event_no")

  selected_records_ls <- map2(
    events,
    event_index,
    function(event, index) {
      if (is.null(event$dataset_name)) {
        data_source <- dataset
      } else {
        data_source <- source_datasets[[event$dataset_name]]
      }
      if (inherits(event, "event")) {
        data_events <- data_source %>%
          group_by(!!!by_vars) %>%
          filter_if(event$condition) %>%
          ungroup()
        if (!is.null(event$mode)) {
          data_events <- filter_extreme(
            data_source,
            by_vars = by_vars,
            order = event$order,
            mode = event$mode
          )
        }
      } else {
        if (is.null(event$order)) {
          event_order <- order
        } else {
          event_order <- event$order
        }
        data_events <- filter_joined(
          dataset,
          by_vars = by_vars,
          join_vars = event$join_vars,
          join_type = event$join_type,
          first_cond = !!event$first_cond,
          order = event_order,
          filter = !!event$condition
        )
      }
      if (is.null(event$keep_vars_source)) {
        event_keep_vars_source <- keep_vars_source
      } else {
        event_keep_vars_source <- event$keep_vars_source
      }
      data_events %>%
        mutate(!!tmp_event_no := index) %>%
        process_set_values_to(set_values_to = event$set_values_to) %>%
        select(!!!event_keep_vars_source, !!!by_vars, names(event$set_values_to))
    }
  )
  selected_records <- bind_rows(selected_records_ls)

  ## tmp obs number within by_vars and a type of event
  tmp_obs <- get_new_tmp_var(selected_records)
  selected_records <- selected_records %>%
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
  new_obs <- selected_records %>%
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
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_extreme_event()`. If the argument is not specified, the input
#'   dataset (`dataset`) of `derive_extreme_event()` is used.
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
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @param set_values_to A named list returned by `exprs()` defining the variables
#'   to be set for the event, e.g. `exprs(PARAMCD = "WSP",
#'   PARAM  = "Worst Sleeping Problems")`. The values can be a symbol, a
#'   character string, a numeric value, `NA` or an expression.
#'
#' @param keep_vars_source Variables to keep from the source dataset
#'
#'   The specified variables are kept for the selected observations. The
#'   variables specified for `by_vars` (of `derive_extreme_event()`) and created
#'   by `set_values_to` are always kept.
#'
#'   *Permitted Values*: A list of expressions where each element is
#'   a symbol or a tidyselect expression, e.g., `exprs(VISIT, VISITNUM,
#'   starts_with("RS"))`.
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_extreme_event()], [event_joined()]
#'
#' @export
#'
#' @return An object of class `event`
event <- function(dataset_name = NULL,
                  condition = NULL,
                  mode = NULL,
                  order = NULL,
                  set_values_to = NULL,
                  keep_vars_source = NULL) {
  out <- list(
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
    keep_vars_source = assert_expr_list(keep_vars_source, optional = TRUE)
  )
  class(out) <- c("event", "event_def", "source", "list")
  out
}

#' Create a `event_joined` Object
#'
#' @description
#'
#' The `event_joined` object is used to define events as input for the
#' `derive_extreme_event()` function. This object should be used if the event
#' does not depend on a single observation of the source dataset but on multiple
#' observations. For example, if the event needs to be confirmed by a second
#' observation of the source dataset.
#'
#' The events are selected by calling `filter_joined()`. See its documentation
#' for more details.
#'
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_extreme_event()`. If the argument is not specified, the input
#'   dataset (`dataset`) of `derive_extreme_event()` is used.
#'
#' @param condition An unquoted condition for selecting the observations, which
#'   will contribute to the extreme event.
#'
#'   *Permitted Values*: an unquoted condition
#'
#' @param join_vars Variables to keep from joined dataset
#'
#'   The variables needed from the other observations should be specified for
#'   this parameter. The specified variables are added to the joined dataset
#'   with suffix ".join". For example to select all observations with `AVALC ==
#'   "Y"` and `AVALC == "Y"` for at least one subsequent visit `join_vars =
#'   exprs(AVALC, AVISITN)` and `filter = AVALC == "Y" & AVALC.join == "Y" &
#'   AVISITN < AVISITN.join` could be specified.
#'
#'   The `*.join` variables are not included in the output dataset.
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
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the subsequent observations, all
#'   observations are removed.
#'
#' @param order If specified, the specified variables or expressions are used to
#'   select the first observation.
#'
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @inheritParams event
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_extreme_event()], [event()]
#'
#' @export
#'
#' @return An object of class `event_joined`
event_joined <- function(dataset_name = NULL,
                         condition,
                         order = NULL,
                         join_vars,
                         join_type,
                         first_cond = NULL,
                         set_values_to = NULL,
                         keep_vars_source = NULL) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name, optional = TRUE),
    condition = assert_filter_cond(enexpr(condition), optional = TRUE),
    order = assert_expr_list(order, optional = TRUE),
    join_vars = assert_vars(join_vars),
    join_type = assert_character_scalar(
      join_type,
      values = c("before", "after", "all"),
      case_sensitive = FALSE
    ),
    first_cond = assert_filter_cond(enexpr(first_cond), optional = TRUE),
    set_values_to = assert_expr_list(
      set_values_to,
      named = TRUE,
      optional = TRUE
    ),
    keep_vars_source = assert_expr_list(keep_vars_source, optional = TRUE)
  )
  class(out) <- c("event_joined", "event_def", "source", "list")
  out
}
