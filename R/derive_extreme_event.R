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
#'   `filter_joined()`. The `condition` field is passed to the `filter` argument.
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
#'       1. The variable specified for `tmp_event_nr_var` is added and set to
#'       the number of the event.
#'       1. Only the variables specified for the `keep_source_vars` field of the
#'       event, and the by variables (`by_vars`) and the variables created by
#'       `set_values_to` are kept.
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
derive_extreme_event <- function(dataset,
                                 by_vars = NULL,
                                 events,
                                 tmp_event_nr_var = NULL,
                                 order,
                                 mode,
                                 source_datasets = NULL,
                                 ignore_event_order = NULL,
                                 check_type = "warning",
                                 set_values_to,
                                 keep_source_vars = exprs(everything())) {
  # Check input parameters
  assert_data_frame(dataset, optional = TRUE)
  assert_vars(by_vars, optional = TRUE)
  assert_list_of(events, "event_def")
  assert_expr_list(order)
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
    deprecate_warn(
      "1.0.0",
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
  assert_varval_list(set_values_to)
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
          data_events <- filter_extreme(
            data_source,
            by_vars = by_vars,
            order = event_order,
            mode = event$mode
          )
        }
      } else {
        data_events <- filter_joined(
          data_source,
          dataset_add = data_source,
          by_vars = by_vars,
          join_vars = event$join_vars,
          join_type = event$join_type,
          first_cond_lower = !!event$first_cond_lower,
          first_cond_upper = !!event$first_cond_upper,
          order = event_order,
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

  ## filter_extreme
  new_obs <- selected_records %>%
    filter_extreme(
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    ) %>%
    mutate(!!!set_values_to) %>%
    select(-!!tmp_event_nr_var)

  # Create output dataset
  bind_rows(dataset, new_obs)
}
#' Add the Worst or Best Observation for Each By Group as New Variables
#'
#' Add the first available record from `events` for each by group as new
#' variables, all variables of the selected observation are kept. It can be used
#' for selecting the extreme observation from a series of user-defined events.
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
#'   `filter_joined()`. The `condition` field is passed to the `filter` argument.
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
#' @param keep_source_vars Variables to keep from the source dataset
#'
#'   For each event the specified variables are kept from the selected
#'   observations. The variables specified for `by_vars` and created by
#'   `set_values_to` are always kept.
#'
#'   *Permitted Values*: A list of expressions where each element is
#'   a symbol or a tidyselect expression, e.g., `exprs(VISIT, VISITNUM,
#'   starts_with("RS"))`.
#'
#' @param new_vars Variables to add
#'
#'   The specified variables from the events are added to the output
#'   dataset. Variables can be renamed by naming the element, i.e., `new_vars =
#'   exprs(<new name> = <old name>)`.
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
#'       `set_values_to` are kept.
#'   1. All selected observations are bound together.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is selected.
#'   1. The variables specified by the `set_values_to` parameter are added to
#'   the selected observations.
#'   1. The variables are added to input dataset.
#'
#'
#' @return The input dataset with the best or worst observation of each by group
#'   added as new variables.
#'
#' @family der_adsl
#' @keywords der_adsl
#'
#' @seealso [event()], [event_joined()], [derive_extreme_event()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#' library(lubridate)
#'
#' adsl <- tribble(
#'   ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
#'   "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
#'   "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""),
#'   "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""),
#'   "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""),
#'   "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
#' )
#'
#' lb <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
#'   "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
#'   "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
#'   "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
#'   "PILOT01",    "LB", "01-1133",    304, "2013-05-01T10:13",
#'   "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
#'   "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
#'   "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
#'   "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
#'   "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
#'   "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
#' ) %>%
#'   mutate(
#'     ADT = convert_dtc_to_dt(LBDTC)
#'   )
#'
#' derive_vars_extreme_event(
#'   adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       dataset_name = "adsl",
#'       condition = !is.na(DTHDT),
#'       set_values_to = exprs(LSTALVDT = DTHDT, DTHFL = "Y")
#'     ),
#'     event(
#'       dataset_name = "lb",
#'       condition = !is.na(ADT),
#'       order = exprs(ADT),
#'       mode = "last",
#'       set_values_to = exprs(LSTALVDT = ADT, DTHFL = "N")
#'     ),
#'     event(
#'       dataset_name = "adsl",
#'       condition = !is.na(TRTEDT),
#'       order = exprs(TRTEDT),
#'       mode = "last",
#'       set_values_to = exprs(LSTALVDT = TRTEDT, DTHFL = "N")
#'     )
#'   ),
#'   source_datasets = list(adsl = adsl, lb = lb),
#'   order = exprs(USUBJID, LSTALVDT, event_nr),
#'   mode = "last",
#'   keep_source_vars = exprs(STUDYID, USUBJID),
#'   new_vars = exprs(LSTALVDT = LSTALVDT, DTHFL = DTHFL)
#' )
#'
#' # Derive DTHCAUS from AE and DS domain data
#' adsl <- tribble(
#'   ~STUDYID,  ~USUBJID,
#'   "STUDY01", "PAT01",
#'   "STUDY01", "PAT02",
#'   "STUDY01", "PAT03"
#' )
#' ae <- tribble(
#'   ~STUDYID,  ~USUBJID, ~AESEQ, ~AEDECOD,       ~AEOUT,  ~AEDTHDTC,
#'   "STUDY01", "PAT01",  12,     "SUDDEN DEATH", "FATAL", "2021-04-04",
#'   "STUDY01", "PAT01",  13,     "CARDIAC ARREST", "FATAL", "2021-04-03",
#' )
#'
#' ds <- tribble(
#'   ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
#'   "STUDY01", "PAT03", 1, "DEATH", "POST STUDY REPORTING OF DEATH", "2022-03-03"
#' )
#'
#' derive_vars_extreme_event(
#'   adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       dataset_name = "ae",
#'       condition = AEOUT == "FATAL",
#'       set_values_to = exprs(DTHCAUS = AEDECOD, DTHDT = convert_dtc_to_dt(AEDTHDTC)),
#'       order = exprs(DTHDT)
#'     ),
#'     event(
#'       dataset_name = "ds",
#'       condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'       set_values_to = exprs(DTHCAUS = DSTERM, DTHDT = convert_dtc_to_dt(DSSTDTC)),
#'       order = exprs(DTHDT)
#'     )
#'   ),
#'   source_datasets = list(ae = ae, ds = ds),
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(USUBJID, DTHDT, event_nr),
#'   mode = "first",
#'   keep_source_vars = exprs(STUDYID, USUBJID),
#'   new_vars = exprs(DTHCAUS = DTHCAUS, DTHDT = DTHDT)
#' )
derive_vars_extreme_event <- function(dataset,
                                      by_vars = NULL,
                                      events,
                                      tmp_event_nr_var = NULL,
                                      order,
                                      mode,
                                      source_datasets = NULL,
                                      check_type = "warning",
                                      new_vars,
                                      keep_source_vars = by_vars) {
  assert_expr_list(new_vars)

  new_obs <- derive_extreme_event(
    dataset = NULL,
    by_vars = by_vars,
    events = events,
    tmp_event_nr_var = event_nr,
    order = order,
    mode = mode,
    source_datasets = source_datasets,
    check_type = check_type,
    set_values_to = new_vars,
    keep_source_vars = keep_source_vars
  )

  derive_vars_merged(
    dataset,
    dataset_add = new_obs,
    by_vars = by_vars
  )
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
#' @seealso [derive_extreme_event()], [event()]
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
    deprecate_warn(
      "1.0.0",
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
