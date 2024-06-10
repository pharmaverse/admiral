#' Add the Worst or Best Observation for Each By Group as New Records
#'
#' Add the first available record from `events` for each by group as new
#' records, all variables of the selected observation are kept. It can be used
#' for selecting the extreme observation from a series of user-defined events.
#' This distinguishes `derive_extreme_event()` from `derive_extreme_records()`,
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
#'   `filter_joined()`. The `condition` field is passed to the `filter_join` argument.
#'
#' @param tmp_event_nr_var Temporary event number variable
#'
#'   The specified variable is added to all source datasets and is set to the
#'   number of the event before selecting the records of the event.
#'
#'   It can be used in `order` to determine which record should be used if
#'   records from more than one event are selected.
#'
#'   The variable is not included in the output dataset.
#'
#' @param order Sort order
#'
#'   If a particular event from `events` has more than one observation, within
#'   the event and by group, the records are ordered by the specified order.
#'
#'   `r roxygen_order_na_handling()`
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If a particular event from `events` has more than one observation,
#'   `"first"`/`"last"` is used to select the first/last record of this type of
#'   event sorting by `order`.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of `event()`
#'   and `event_joined()` refers to the dataset provided in the list.
#'
#' @param ignore_event_order Ignore event order
#'
#'   `r lifecycle::badge("deprecated")`
#'
#'   This argument is *deprecated*. If event order should be ignored, please
#'   specify neither `ignore_event_order` nor `tmp_event_nr_var`. If the event
#'   order should be considered, don't specify `ignore_event_order` but specify
#'   `tmp_event_nr_var` and and the specified variable to `order`.
#'
#'   If the argument is set to `TRUE`, all events defined by `events` are
#'   considered equivalent. If there is more than one observation per by group
#'   the first or last (with respect to `mode` and `order`) is select without
#'   taking the order of the events into account.
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param keep_source_vars Variables to keep from the source dataset
#'
#'  For each event the specified variables are kept from the selected
#'  observations. The variables specified for `by_vars` and created by
#'  `set_values_to` are always kept. The `keep_source_vars` field of
#'  the event will take precedence over the value of the `keep_source_vars`
#'  argument.
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
#'       1. The variable specified for `tmp_event_nr_var` is added and set to
#'       the number of the event.
#'       1. Only the variables specified for the `keep_source_vars` field of the
#'       event, and the by variables (`by_vars`) and the variables created by
#'       `set_values_to` are kept. If `keep_source_vars = NULL` is used for an event
#'       in `derive_extreme_event()` the value of the `keep_source_vars` argument of
#'       `derive_extreme_event()` is used.
#'   1. All selected observations are bound together.
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
#' @seealso [event()], [event_joined()], [derive_vars_extreme_event()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#' library(lubridate)
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
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(event_nr, desc(ADY)),
#'   mode = "first",
#'   set_values_to = exprs(
#'     PARAMCD = "WSP",
#'     PARAM = "Worst Sleeping Problems"
#'   )
#' )
#'
#' # Use different mode by event
#' adhy <- tribble(
#'   ~USUBJID, ~AVISITN, ~CRIT1FL,
#'   "1",             1, "Y",
#'   "1",             2, "Y",
#'   "2",             1, "Y",
#'   "2",             2, NA_character_,
#'   "2",             3, "Y",
#'   "2",             4, NA_character_
#' ) %>%
#'   mutate(
#'     PARAMCD = "ALKPH",
#'     PARAM = "Alkaline Phosphatase (U/L)"
#'   )
#'
#' derive_extreme_event(
#'   adhy,
#'   by_vars = exprs(USUBJID),
#'   events = list(
#'     event(
#'       condition = is.na(CRIT1FL),
#'       set_values_to = exprs(AVALC = "N")
#'     ),
#'     event(
#'       condition = CRIT1FL == "Y",
#'       mode = "last",
#'       set_values_to = exprs(AVALC = "Y")
#'     )
#'   ),
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(event_nr, AVISITN),
#'   mode = "first",
#'   keep_source_vars = exprs(AVISITN),
#'   set_values_to = exprs(
#'     PARAMCD = "ALK2",
#'     PARAM = "ALKPH <= 2 times ULN"
#'   )
#' )
#'
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
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(event_nr, ADT),
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
#'
derive_extreme_event <- function(dataset = NULL,
                                 by_vars,
                                 events,
                                 tmp_event_nr_var = NULL,
                                 order,
                                 mode,
                                 source_datasets = NULL,
                                 ignore_event_order = NULL,
                                 check_type = "warning",
                                 set_values_to = NULL,
                                 keep_source_vars = exprs(everything())) {
  # Check input parameters
  assert_data_frame(dataset, optional = TRUE)
  assert_vars(by_vars)
  assert_list_of(events, "event_def")
  assert_expr_list(order)
  mode <- assert_character_scalar(mode, values = c("first", "last"), case_sensitive = FALSE)
  assert_list_of(source_datasets, "data.frame")
  source_names <- names(source_datasets)
  assert_list_element(
    list = events,
    element = "dataset_name",
    condition = map_lgl(dataset_name, is.null) | dataset_name %in% source_names,
    source_names = source_names,
    message_text = c(
      paste0(
        "The dataset names must be included in the list specified for the ",
        "{.arg source_datasets} argument."
      ),
      i = paste(
        "Following names were provided by {.arg source_datasets}:",
        ansi_collapse(source_names)
      )
    )
  )

  if (!is.null(ignore_event_order)) {
    if (ignore_event_order) {
      deprecate_details <- paste(
        "The event order is ignored by default.",
        "Please remove `ignore_event_order = TRUE` from the call.",
        sep = "\n"
      )
    } else {
      deprecate_details <- c(
        "Please remove `ignore_event_order = FALSE` from the call.",
        "Specify `tmp_event_nr_var`.",
        "Add the specified variable to `order`."
      )
    }
    deprecate_stop(
      "1.1.0",
      "derive_extreme_event(ignore_event_order=)",
      "derive_extreme_event(tmp_event_nr_var=)",
      details = deprecate_details
    )
  }
  tmp_event_nr_var <- assert_symbol(enexpr(tmp_event_nr_var), optional = TRUE)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_varval_list(set_values_to, optional = TRUE)
  keep_source_vars <- assert_expr_list(keep_source_vars)

  # Create new observations
  ## Create a dataset (selected_records) from `events`
  event_index <- as.list(seq_along(events))

  selected_records_ls <- map2(
    events,
    event_index,
    function(event, index) {
      if (is.null(event$dataset_name)) {
        data_source <- dataset
      } else {
        data_source <- source_datasets[[event$dataset_name]]
      }
      if (!is.null(tmp_event_nr_var)) {
        data_source <- mutate(data_source, !!tmp_event_nr_var := index)
      }
      if (is.null(event$order)) {
        event_order <- order
      } else {
        event_order <- event$order
      }
      if (inherits(event, "event")) {
        data_events <- data_source %>%
          group_by(!!!by_vars) %>%
          filter_if(event$condition) %>%
          ungroup()
        if (!is.null(event$mode)) {
          if (check_type != "none") {
            # Check for duplicates
            signal_duplicate_records(
              dataset = data_events,
              by_vars = append(by_vars, event_order),
              msg = paste(
                "Check duplicates: ", event$dataset_name,
                "dataset contains duplicate records with respect to",
                "{.var {replace_values_by_names(by_vars)}}"
              ),
              cnd_type = check_type
            )
          }

          data_events <- filter_extreme(
            data_events,
            by_vars = by_vars,
            order = event_order,
            mode = event$mode,
            check_type = "none"
          )
        }
      } else {
        if (check_type != "none") {
          # Check for duplicates
          signal_duplicate_records(
            dataset = data_source,
            by_vars = append(by_vars, event_order),
            msg = paste(
              "Check duplicates: ", event$dataset_name,
              "dataset contains duplicate records with respect to",
              "{.var {replace_values_by_names(by_vars)}}"
            ),
            cnd_type = check_type
          )
        }

        data_events <- filter_joined(
          data_source,
          dataset_add = data_source,
          by_vars = by_vars,
          join_vars = event$join_vars,
          join_type = event$join_type,
          first_cond_lower = !!event$first_cond_lower,
          first_cond_upper = !!event$first_cond_upper,
          order = event_order,
          check_type = "none",
          filter_join = !!event$condition
        )
      }
      if (is.null(event$keep_source_vars)) {
        event_keep_source_vars <- keep_source_vars
      } else {
        event_keep_source_vars <- event$keep_source_vars
      }
      data_events %>%
        process_set_values_to(set_values_to = event$set_values_to) %>%
        select(
          !!!event_keep_source_vars, !!tmp_event_nr_var, !!!by_vars,
          names(event$set_values_to)
        )
    }
  )
  selected_records <- bind_rows(selected_records_ls)

  if (check_type != "none") {
    # Check for duplicates
    signal_duplicate_records(
      dataset = selected_records,
      by_vars = append(by_vars, order),
      msg = paste(
        "Check duplicates: the dataset which consists of all records selected",
        "for any of the events defined by {.arg events} contains duplicate records",
        "with respect to {.var {replace_values_by_names(by_vars)}}"
      ),
      cnd_type = check_type
    )
  }

  ## filter_extreme
  new_obs <- selected_records %>%
    filter_extreme(
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = "none"
    ) %>%
    mutate(!!!set_values_to) %>%
    select(-!!tmp_event_nr_var)

  # Create output dataset
  bind_rows(dataset, new_obs)
}
