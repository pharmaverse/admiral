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
#' @permitted [event]
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
#' @permitted [var]
#'
#' @param order Sort order
#'
#'   If a particular event from `events` has more than one observation, within
#'   the event and by group, the records are ordered by the specified order.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [var_list]
#'
#' @param mode Selection mode (first or last)
#'
#'   If a particular event from `events` has more than one observation,
#'   `"first"`/`"last"` is used to select the first/last record of this type of
#'   event sorting by `order`.
#'
#' @permitted [mode]
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of `event()`
#'   and `event_joined()` refers to the dataset provided in the list.
#'
#' @permitted [dataset_list]
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations.
#'
#'   Set a list of variables to some specified value for the new records
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value, an expression or NA.
#'
#'   For example:
#'   ```
#'     set_values_to = exprs(
#'       PARAMCD = "WOBS",
#'       PARAM = "Worst Observations"
#'     )
#'   ```
#'
#' @permitted [expr_list_formula]
#'
#' @param keep_source_vars Variables to keep from the source dataset
#'
#'  For each event the specified variables are kept from the selected
#'  observations. The variables specified for `by_vars` and created by
#'  `set_values_to` are always kept. The `keep_source_vars` field of
#'  the event will take precedence over the value of the `keep_source_vars`
#'  argument.
#'
#' @permitted [var_list]
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
#' `r roxygen_save_memory()`
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
#' @examplesx
#'
#' @caption Add a new record for the worst observation using `event()` objects
#' @info For each subject, the observation containing the worst sleeping problem
#'   (if any exist) should be identified and added as a new record, retaining
#'   all variables from the original observation. If multiple occurrences of the
#'   worst sleeping problem occur, or no sleeping problems, then take the
#'   observation occurring at the latest day.
#'
#' - The groups for which new records are added are specified by the `by_vars`
#'   argument. Here for each *subject* a record should be added. Thus
#'   `by_vars = exprs(STUDYID, USUBJID)` is specified.
#' - The sets of possible sleeping problems are passed through the `events`
#'   argument as `event()` objects. Each event contains a `condition` which
#'   may or may not be satisfied by each record (or possibly a group of
#'   records) within the input dataset `dataset`. Summary functions such as
#'   `any()` and `all()` are often handy to use within conditions, as is done
#'   here for the third event, which checks that the subject had no sleeping
#'   issues. The final event uses a catch-all `condition = TRUE` to ensure all
#'   subjects have a new record derived. Note that in this example, as no
#'   condition involves analysis of __cross-comparison values of within  records__,
#'   it is sufficient to use `event()` objects rather than `event_joined()` -
#'   see the next example for a more complex condition.
#' - If any subject has one or more records satisfying the conditions from
#'   events, we can select just one record using the `order` argument. In this
#'   example, the first argument passed to `order` is `event_nr`, which is a
#'   temporary variable created through the `tmp_event_nr_var` argument, which
#'   numbers the events consecutively. Since `mode = "first"`, we only consider
#'   the first event for which a condition is satisfied. Within that event, we
#'   consider only the observation with the latest day, because the second
#'   argument for the order is `desc(ADY)`.
#' - Once a record is identified as satisfying an event's condition, a new
#'   observation is created by the following process:
#'   \enumerate{
#'     \item the selected record is copied,
#'     \item the variables specified in the event's `set_values_to` (here,
#'     `AVAL` and `AVALC`) are created/updated,
#'     \item the variables specified in `keep_source_vars` (here, `ADY` does due
#'     to the use of the tidyselect expression `everything()`) (plus `by_vars`
#'     and the variables from `set_values_to`) are kept,
#'     \item the variables specified in the global `set_values_to` (here,
#'     `PARAM` and `PARAMCD`) are created/updated.
#'   }
#'
#' @code
#' library(tibble, warn.conflicts = FALSE)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' adqs1 <- tribble(
#'   ~USUBJID, ~PARAMCD,         ~AVALC,        ~ADY,
#'   "1",      "NO SLEEP",       "N",              1,
#'   "1",      "WAKE UP 3X",     "N",              2,
#'   "2",      "NO SLEEP",       "N",              1,
#'   "2",      "WAKE UP 3X",     "Y",              2,
#'   "2",      "WAKE UP 3X",     "Y",              3,
#'   "3",      "NO SLEEP",       NA_character_,    1
#' ) %>%
#' mutate(STUDYID = "AB42")
#'
#' derive_extreme_event(
#'   adqs1,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
#'     ),
#'     event(
#'       condition = PARAMCD == "WAKE UP 3X" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "Waking up three times", AVAL = 2)
#'     ),
#'     event(
#'       condition = all(AVALC == "N"),
#'       set_values_to = exprs(AVALC = "No sleeping problems", AVAL = 3)
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
#'     PARAM = "Worst Sleeping Problem"
#'   ),
#'   keep_source_vars = exprs(everything())
#' ) %>%
#' select(-STUDYID)
#'
#' @caption Events based on comparison across records (`event_joined()`)
#' @info We'll now extend the above example. Specifically, we consider a new
#'    possible worst sleeping problem, namely if a subject experiences no
#'    sleep on consecutive days.
#'
#'  - The "consecutive days" portion of the condition requires records to be
#'    compared with each other. This is done by using an `event_joined()` object,
#'    specifically by passing `dataset_name = adqs2` to it so that the `adqs2`
#'    dataset is joined onto itself. The `condition` now checks for two
#'    no sleep records, and crucially compares the `ADY` values to see if
#'    they differ by one day. The `.join` syntax distinguishes between the
#'    `ADY` value of the parent and joined datasets. As the condition involves
#'    `AVALC`, `PARAMCD` and `ADY`, we specify these variables with `join_vars`,
#'    and finally, because we wish to compare all records with each other, we
#'    select `join_type = "all"`.
#'
#' @code
#' adqs2 <- tribble(
#'    ~USUBJID, ~PARAMCD,     ~AVALC, ~ADY,
#'    "4",      "WAKE UP",    "N",    1,
#'    "4",      "NO SLEEP",   "Y",    2,
#'    "4",      "NO SLEEP",   "Y",    3,
#'    "5",      "NO SLEEP",   "N",    1,
#'    "5",      "NO SLEEP",   "Y",    2,
#'    "5",      "WAKE UP 3X", "Y",    3,
#'    "5",      "NO SLEEP",   "Y",    4
#' ) %>%
#' mutate(STUDYID = "AB42")
#'
#' derive_extreme_event(
#'   adqs2,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event_joined(
#'       join_vars = exprs(AVALC, PARAMCD, ADY),
#'       join_type = "all",
#'       condition = PARAMCD == "NO SLEEP" & AVALC == "Y" &
#'         PARAMCD.join == "NO SLEEP" & AVALC.join == "Y" &
#'         ADY == ADY.join + 1,
#'       set_values_to = exprs(AVALC = "No sleep two nights in a row", AVAL = 0)
#'     ),
#'     event(
#'       condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
#'     ),
#'     event(
#'       condition = PARAMCD == "WAKE UP 3X" & AVALC == "Y",
#'       set_values_to = exprs(AVALC = "Waking up three times", AVAL = 2)
#'     ),
#'     event(
#'       condition = all(AVALC == "N"),
#'       set_values_to = exprs(
#'         AVALC = "No sleeping problems", AVAL = 3
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
#'     PARAM = "Worst Sleeping Problem"
#'   ),
#'   keep_source_vars = exprs(everything())
#' ) %>%
#' select(-STUDYID)
#'
#' @caption Specifying different arguments across `event()` objects
#' @info Here we consider a Hy's Law use case. We are interested in
#'   knowing whether a subject's Alkaline Phosphatase has ever been
#'   above twice the upper limit of normal range. If so, i.e. if
#'   `CRIT1FL` is `Y`, we are interested in the record for the first
#'   time this occurs, and if not, we wish to retain the last record.
#'   As such, for this case now we need to vary our usage of the
#'   `mode` argument dependent on the `event()`.
#'
#'  - In first `event()`, since we simply seek the first time that
#'    `CRIT1FL` is `"Y"`, it's enough to specify the `condition`,
#'    because we inherit `order` and `mode` from the main
#'    `derive_extreme_event()` call here which will automatically
#'    select the first occurrence by `AVISITN`.
#'  - In the second `event()`, we select the last record among the
#'    full set of records where `CRIT1FL` are all `"N"` by additionally
#'    specifying `mode = "last"` within the `event()`.
#'  - Note now the usage of `keep_source_vars = exprs(AVISITN)`
#'    rather than `everything()` as in the previous example. This
#'    is done to ensure `CRIT1` and `CRIT1FL` are not populated for
#'    the new records.
#'
#' @code
#' adhy <- tribble(
#'   ~USUBJID, ~AVISITN,              ~CRIT1, ~CRIT1FL,
#'   "1",             1, "ALT > 2 times ULN", "N",
#'   "1",             2, "ALT > 2 times ULN", "N",
#'   "2",             1, "ALT > 2 times ULN", "N",
#'   "2",             2, "ALT > 2 times ULN", "Y",
#'   "2",             3, "ALT > 2 times ULN", "N",
#'   "2",             4, "ALT > 2 times ULN", "Y"
#' ) %>%
#'   mutate(
#'     PARAMCD = "ALT",
#'     PARAM = "ALT (U/L)",
#'     STUDYID = "AB42"
#'   )
#'
#' derive_extreme_event(
#'   adhy,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       condition = CRIT1FL == "Y",
#'       set_values_to = exprs(AVALC = "Y")
#'     ),
#'     event(
#'       condition = CRIT1FL == "N",
#'       mode = "last",
#'       set_values_to = exprs(AVALC = "N")
#'     )
#'   ),
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(event_nr, AVISITN),
#'   mode = "first",
#'   keep_source_vars = exprs(AVISITN),
#'   set_values_to = exprs(
#'     PARAMCD = "ALT2",
#'     PARAM = "ALT > 2 times ULN"
#'   )
#' ) %>%
#'   select(-STUDYID)
#'
#' @caption A more complex example: Confirmed Best Overall Response
#' (`first/last_cond_upper`, `join_type`, `source_datasets`)
#' @info The final example showcases a use of `derive_extreme_event()`
#'   to calculate the Confirmed Best Overall Response (CBOR) in an
#'   `ADRS` dataset, as is common in many oncology trials. This example
#'   builds on all the previous ones and thus assumes a baseline level
#'   of confidence with `derive_extreme_event()`.
#'
#'   The following `ADSL` and `ADRS` datasets will be used
#'   throughout:
#'
#' @code
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
#' mutate(
#'   TRTSDT = ymd(TRTSDTC),
#'   STUDYID = "AB42"
#' )
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
#'     STUDYID = "AB42",
#'     PARAMCD = "OVR",
#'     PARAM = "Overall Response by Investigator"
#'   ) %>%
#'   derive_vars_merged(
#'     dataset_add = adsl,
#'     by_vars = exprs(STUDYID, USUBJID),
#'     new_vars = exprs(TRTSDT)
#'   )
#'
#' @info Since the CBOR derivation contains multiple complex parts, it's
#'   convenient to make use of the `description` argument within each event object
#'   to describe what condition is being checked.
#'
#'   - For the Confirmed Response (CR), for each `"CR"` record in the original `ADRS`
#'     dataset that will be identified by the first part of the `condition` argument
#'     (`AVALC == "CR"`), we need to use the `first_cond_upper` argument to limit the
#'     group of observations to consider alongside it. Namely, we need to look up to
#'     and including the second CR (`AVALC.join == "CR"`) over 28 days from the first
#'     one (`ADT.join >= ADT + 28`). The observations satisfying `first_cond_upper`
#'     then form part of our "join group", meaning that the remaining portions of
#'     `condition` which reference joined variables are limited to this group.
#'     In particular, within `condition` we use `all()` to check that all observations
#'     are either `"CR"` or `"NE"`, and `count_vals()` to ensure at most one is
#'     `"NE"`.
#'
#'     Note that the selection of `join_type = "after"` is critical here, due to the
#'     fact that the restriction implied by `join_type` is applied before the one
#'     implied by `first_cond_upper`. Picking the first subject (who was correctly
#'     identified as a confirmed responder) as an example, selecting
#'     `join_type = "all"` instead of `"after"` would mean the first `"PR"` record
#'     from `"2020-01-01"` would also be considered when evaluating the
#'     `all(AVALC.join %in% c("CR", "NE"))` portion of `condition`. In turn, the
#'     condition would not be satisfied anymore, and in this case, following the
#'     later event logic shows the subject would be considered a partial responder
#'     instead.
#'
#'  - The Partial Response (PR), is very similar; with the difference being that the
#'    first portion of `condition` now references `"PR"` and `first_cond_upper`
#'    accepts a confirmatory `"PR"` or `"CR"` 28 days later. Note that now we must add
#'    `"PR"` as an option within the `all()` condition to account for confirmatory
#'    `"PR"`s.
#'
#'  - The Stable Disease (SD), Progressive Disease (PD) and Not Evaluable (NE)
#'    events are simpler and just require `event()` calls.
#'
#'  - Finally, we use a catch-all `event()` with `condition = TRUE` and
#'    `dataset_name = "adsl"` to identify those subjects who do not appear in `ADRS`
#'    and list their CBOR as `"MISSING"`. Note here the fact that `dataset_name` is
#'    set to `"adsl"`, which is a new source dataset. As such it's important in the
#'    main `derive_extreme_event()` call to list `adsl` as another source dataset
#'    with `source_datasets = list(adsl = adsl)`.
#'
#' @code
#' derive_extreme_event(
#'   adrs,
#'   by_vars = exprs(STUDYID, USUBJID),
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
#'       first_cond_upper = AVALC.join == "CR" & ADT.join >= ADT + 28,
#'       condition = AVALC == "CR" &
#'         all(AVALC.join %in% c("CR", "NE")) &
#'         count_vals(var = AVALC.join, val = "NE") <= 1,
#'       set_values_to = exprs(AVALC = "CR")
#'     ),
#'     event_joined(
#'       description = paste(
#'         "PR needs to be confirmed by a second CR or PR at least 28 days later,",
#'         "at most one NE is acceptable between the two assessments"
#'       ),
#'       join_vars = exprs(AVALC, ADT),
#'       join_type = "after",
#'       first_cond_upper = AVALC.join %in% c("CR", "PR") & ADT.join >= ADT + 28,
#'       condition = AVALC == "PR" &
#'         all(AVALC.join %in% c("CR", "PR", "NE")) &
#'         count_vals(var = AVALC.join, val = "NE") <= 1,
#'       set_values_to = exprs(AVALC = "PR")
#'     ),
#'     event(
#'       description = paste(
#'         "CR, PR, or SD are considered as SD if occurring at least 28",
#'         "after treatment start"
#'       ),
#'       condition = AVALC %in% c("CR", "PR", "SD") & ADT >= TRTSDT + 28,
#'       set_values_to = exprs(AVALC = "SD")
#'     ),
#'     event(
#'       condition = AVALC == "PD",
#'       set_values_to = exprs(AVALC = "PD")
#'     ),
#'     event(
#'       condition = AVALC %in% c("CR", "PR", "SD", "NE"),
#'       set_values_to = exprs(AVALC = "NE")
#'     ),
#'     event(
#'       description = "Set response to MISSING for patients without records in ADRS",
#'       dataset_name = "adsl",
#'       condition = TRUE,
#'       set_values_to = exprs(AVALC = "MISSING"),
#'       keep_source_vars = exprs(TRTSDT)
#'     )
#'   ),
#'   set_values_to = exprs(
#'     PARAMCD = "CBOR",
#'     PARAM = "Best Confirmed Overall Response by Investigator"
#'   )
#' ) %>%
#'   filter(PARAMCD == "CBOR") %>%
#'   select(-STUDYID, -ADTC)
#'
#' @caption Further examples
#' @info Equivalent examples for using the`check_type` argument can be found in
#'   `derive_extreme_records()`.
derive_extreme_event <- function(dataset = NULL,
                                 by_vars,
                                 events,
                                 tmp_event_nr_var = NULL,
                                 order,
                                 mode,
                                 source_datasets = NULL,
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
