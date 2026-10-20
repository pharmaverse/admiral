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
#' @permitted list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If a particular event from `events` has more than one observation,
#'   `"first"`/`"last"` is used to select the first/last record of this type of
#'   event sorting by `order`.
#'
#' @permitted `"first"`, `"last"`
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of `event()`
#'   and `event_joined()` refers to the dataset provided in the list.
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
#'   1. All selected observations are bound together.
#'   1. For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is selected.
#'   1. The variables specified by the `new_vars` parameter are added to
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
#' @examplesx
#'
#' @caption Derive the last known alive date from multiple sources
#'
#' @info All records in these examples are synthetic. For each subject, select
#'   the latest of the death date, a dated laboratory assessment, and treatment
#'   end. These sources are an illustrative study rule; their eligibility should
#'   be specified for the study. Keep the selected date and its source in ADSL.
#'
#' - `dataset` is the ADSL dataset to receive the new variables. `by_vars`
#'   identifies subjects by study and subject ID, so records from different
#'   studies are not combined even if they share a subject ID.
#' - `source_datasets` supplies the datasets named by each `event()`.
#' - The laboratory event first selects the last nonmissing `ADT` per subject
#'   using its own `order` and `mode`. The outer `order` and `mode` then select
#'   the latest date across all events.
#' - `tmp_event_nr_var` numbers events in their listed order. With
#'   `order = exprs(LSTALVDT, desc(event_nr))` and `mode = "last"`, an equal-date
#'   tie favors the event listed first: death, then laboratory, then treatment
#'   end. Subject `09-1088` illustrates a death/laboratory tie.
#' - `new_vars` selects the variables to merge into ADSL. All ADSL subjects are
#'   retained; subject `09-1099` has no eligible source record, so both new
#'   variables are missing. The temporary event number is not retained.
#' - `check_type = "error"` stops if candidate records are not unique with
#'   respect to the grouping and ordering variables. The event number resolves
#'   equal-date ties across sources in this example; `check_type` detects
#'   duplicates but does not resolve them.
#'
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' adsl <- tribble(
#'   ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
#'   "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
#'   "PILOT01", "01-1133", ymd("2013-04-28"), ymd(NA_character_),
#'   "PILOT01", "01-1211", ymd("2013-01-12"), ymd(NA_character_),
#'   "PILOT01", "09-1081", ymd("2014-04-27"), ymd(NA_character_),
#'   "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
#'   "PILOT01", "09-1099", ymd(NA_character_), ymd(NA_character_)
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
#'   "PILOT01",    "LB", "09-1088",    322, "2014-11-01T13:25"
#' ) %>%
#'   mutate(
#'     ADT = convert_dtc_to_dt(LBDTC)
#'   )
#'
#' derive_vars_extreme_event(
#'   dataset = adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       dataset_name = "adsl",
#'       condition = !is.na(DTHDT),
#'       set_values_to = exprs(LSTALVDT = DTHDT, LALVSRC = "Death date")
#'     ),
#'     event(
#'       dataset_name = "lb",
#'       condition = !is.na(ADT),
#'       order = exprs(ADT),
#'       mode = "last",
#'       set_values_to = exprs(LSTALVDT = ADT, LALVSRC = "Laboratory")
#'     ),
#'     event(
#'       dataset_name = "adsl",
#'       condition = !is.na(TRTEDT),
#'       order = exprs(TRTEDT),
#'       mode = "last",
#'       set_values_to = exprs(LSTALVDT = TRTEDT, LALVSRC = "Treatment end")
#'     )
#'   ),
#'   source_datasets = list(adsl = adsl, lb = lb),
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(LSTALVDT, desc(event_nr)),
#'   mode = "last",
#'   check_type = "error",
#'   new_vars = exprs(LSTALVDT, LALVSRC)
#' )
#'
#' @caption Derive a reported death cause with an explicit tie-break rule
#'
#' @info Select the earliest dated qualifying report from AE or DS and add its
#'   date and reported cause to ADSL. This is an illustrative source-selection
#'   rule, not clinical adjudication of conflicting reports.
#'
#' - Each event maps its source variables to the common names `DTHCAUS` and
#'   `DTHDT`. The final `order` uses these names across the combined events.
#' - Date takes precedence over source: `order = exprs(DTHDT, event_nr)` and
#'   `mode = "first"` choose the earliest date, with AE preferred over DS only
#'   when dates are equal. Subject `PAT04` illustrates this tie.
#' - `check_type = "error"` stops if multiple eligible reports from the same
#'   source share a subject and date. Resolve such ambiguity using a
#'   study-specific rule before deriving the reported cause.
#' - `DTHCAUS` and `DTHDT` are taken together from the selected qualifying
#'   record. `DTHDT` is not an independent derivation from all death records.
#'   The DS record for `PAT03` reports death on 2022-03-03, but its text does
#'   not satisfy this example's cause-selection condition, so both derived
#'   variables are missing. This does not mean the subject is alive or that
#'   their death date is unknown. Derive the death date separately if records
#'   without an eligible reported cause must also contribute to it.
#'
#' @code
#' adsl <- tribble(
#'   ~STUDYID,  ~USUBJID,
#'   "STUDY01", "PAT01",
#'   "STUDY01", "PAT02",
#'   "STUDY01", "PAT03",
#'   "STUDY01", "PAT04"
#' )
#' ae <- tribble(
#'   ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
#'   "STUDY01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04",
#'   "STUDY01", "PAT01", 13, "CARDIAC ARREST", "FATAL", "2021-04-03",
#'   "STUDY01", "PAT04",  1, "RESPIRATORY FAILURE", "FATAL", "2022-02-01"
#' )
#'
#' ds <- tribble(
#'   ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED",
#'     "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", 2, "RANDOMIZATION",
#'     "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", 3, "DEATH",
#'     "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
#'   "STUDY01", "PAT03", 1, "DEATH",
#'     "POST STUDY REPORTING OF DEATH", "2022-03-03",
#'   "STUDY01", "PAT04", 1, "DEATH",
#'     "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01"
#' )
#'
#' derive_vars_extreme_event(
#'   adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event(
#'       dataset_name = "ae",
#'       condition = AEOUT == "FATAL",
#'       set_values_to = exprs(
#'         DTHCAUS = AEDECOD,
#'         DTHDT = convert_dtc_to_dt(AEDTHDTC)
#'       )
#'     ),
#'     event(
#'       dataset_name = "ds",
#'       condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
#'       set_values_to = exprs(
#'         DTHCAUS = DSTERM,
#'         DTHDT = convert_dtc_to_dt(DSSTDTC)
#'       )
#'     )
#'   ),
#'   source_datasets = list(ae = ae, ds = ds),
#'   tmp_event_nr_var = event_nr,
#'   order = exprs(DTHDT, event_nr),
#'   mode = "first",
#'   check_type = "error",
#'   new_vars = exprs(DTHCAUS, DTHDT)
#' )
#'
#' @caption Derive a confirmation date using `event_joined()`
#'
#' @info For this synthetic study, a laboratory elevation is confirmed when ALT
#'   exceeds three times the upper limit of normal (ULN) on two consecutive
#'   calendar days. Derive the first confirmation date, i.e., the second date
#'   of the pair. This illustrates a study-specific programming rule, not a
#'   diagnosis. `AVAL` is the ALT result in U/L and `ANRHI` is its ULN in U/L.
#'
#' - The condition compares each assessment with earlier assessments from the
#'   same subject using `event_joined()`, `join_type = "before"`, and
#'   `order = exprs(ADT)`. Variables from earlier records have the `.join`
#'   suffix. `ADT == ADT.join + 1` requires consecutive dates, not merely
#'   consecutive records. The example has one ALT assessment per subject/day.
#' - Both results must exceed their respective ULN threshold. Missing values
#'   do not establish confirmation. The selected record supplies `ADT` from the
#'   confirming assessment; the outer `mode = "first"` selects the first such date.
#' - `new_vars` renames `ADT` to `ALTCDT`. Subjects with a gap between elevated
#'   results (`PAT02`), an unconfirmed elevation (`PAT04`), or no laboratory
#'   records (`PAT05`) remain in ADSL with a missing confirmation date.
#'
#' @code
#' adsl <- tribble(
#'   ~STUDYID,  ~USUBJID,
#'   "STUDY01", "PAT01",
#'   "STUDY01", "PAT02",
#'   "STUDY01", "PAT03",
#'   "STUDY01", "PAT04",
#'   "STUDY01", "PAT05"
#' )
#'
#' adlb <- tribble(
#'   ~STUDYID,  ~USUBJID, ~ADT,         ~AVAL, ~ANRHI,
#'   "STUDY01", "PAT01",  "2024-01-01",  160,     40,
#'   "STUDY01", "PAT01",  "2024-01-02",  200,     40,
#'   "STUDY01", "PAT01",  "2024-01-03",  240,     40,
#'   "STUDY01", "PAT02",  "2024-01-01",  160,     40,
#'   "STUDY01", "PAT02",  "2024-01-03",  200,     40,
#'   "STUDY01", "PAT03",  "2024-01-01",   80,     40,
#'   "STUDY01", "PAT03",  "2024-01-02",  160,     40,
#'   "STUDY01", "PAT03",  "2024-01-03",  200,     40,
#'   "STUDY01", "PAT04",  "2024-01-01",  160,     40,
#'   "STUDY01", "PAT04",  "2024-01-02",   NA,     40
#' ) %>%
#'   mutate(ADT = ymd(ADT))
#'
#' derive_vars_extreme_event(
#'   adsl,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   events = list(
#'     event_joined(
#'       dataset_name = "adlb",
#'       order = exprs(ADT),
#'       join_vars = exprs(ADT, AVAL, ANRHI),
#'       join_type = "before",
#'       condition = AVAL > 3 * ANRHI & AVAL.join > 3 * ANRHI.join &
#'         ADT == ADT.join + 1
#'     )
#'   ),
#'   source_datasets = list(adlb = adlb),
#'   order = exprs(ADT),
#'   mode = "first",
#'   new_vars = exprs(ALTCDT = ADT)
#' )
#'
#' @caption Further examples
#'
#' @info See the "Handling duplicates" section of [derive_extreme_records()]
#'   for examples of the default `check_type = "warning"`, inspecting duplicate
#'   records with [get_duplicates_dataset()], and using `check_type = "none"`.
#'   Disabling the check does not resolve ties; an explicit selection rule is
#'   needed if tied records would give different derived values.
#'
#'   See [derive_extreme_event()] for further examples of `event()` and
#'   `event_joined()`, including event-specific grouping, ordering, and source
#'   restrictions. Here those events supply variables to `dataset` rather than
#'   adding records.
derive_vars_extreme_event <- function(dataset,
                                      by_vars,
                                      events,
                                      tmp_event_nr_var = NULL,
                                      order,
                                      mode,
                                      source_datasets = NULL,
                                      check_type = "warning",
                                      new_vars) {
  assert_expr_list(new_vars)
  tmp_event_nr_var <- assert_symbol(enexpr(tmp_event_nr_var), optional = TRUE)

  new_obs <- derive_extreme_event(
    dataset = NULL,
    by_vars = by_vars,
    events = events,
    tmp_event_nr_var = !!tmp_event_nr_var,
    order = order,
    mode = mode,
    source_datasets = source_datasets,
    check_type = check_type
  )

  derive_vars_merged(
    dataset,
    dataset_add = new_obs,
    new_vars = new_vars,
    by_vars = by_vars
  ) %>%
    as_admiral_df()
}
