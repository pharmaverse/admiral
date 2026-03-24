# Add the Worst or Best Observation for Each By Group as New Records

Add the first available record from `events` for each by group as new
records, all variables of the selected observation are kept. It can be
used for selecting the extreme observation from a series of user-defined
events. This distinguishes `derive_extreme_event()` from
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
where extreme records are derived based on certain order of existing
variables.

## Usage

``` r
derive_extreme_event(
  dataset = NULL,
  by_vars,
  events,
  tmp_event_nr_var = NULL,
  order,
  mode,
  source_datasets = NULL,
  check_type = "warning",
  set_values_to = NULL,
  keep_source_vars = exprs(everything())
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- events:

  Conditions and new values defining events

  A list of
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  or
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  objects is expected. Only observations listed in the `events` are
  considered for deriving extreme event. If multiple records meet the
  filter `condition`, take the first record sorted by `order`. The data
  is grouped by `by_vars`, i.e., summary functions like
  [`all()`](https://rdrr.io/r/base/all.html) or
  [`any()`](https://rdrr.io/r/base/any.html) can be used in `condition`.

  For
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  events the observations are selected by calling
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md).
  The `condition` field is passed to the `filter_join` argument.

  Permitted values

  :   an
      [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
      or
      [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
      object

  Default value

  :   none

- tmp_event_nr_var:

  Temporary event number variable

  The specified variable is added to all source datasets and is set to
  the number of the event before selecting the records of the event.

  It can be used in `order` to determine which record should be used if
  records from more than one event are selected.

  The variable is not included in the output dataset.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- order:

  Sort order

  If a particular event from `events` has more than one observation,
  within the event and by group, the records are ordered by the
  specified order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If a particular event from `events` has more than one observation,
  `"first"`/`"last"` is used to select the first/last record of this
  type of event sorting by `order`.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- source_datasets:

  Source datasets

  A named list of datasets is expected. The `dataset_name` field of
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  and
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  refers to the dataset provided in the list.

  Permitted values

  :   named list of datasets, e.g., `list(adsl = adsl, ae = ae)`

  Default value

  :   `NULL`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the input dataset are not unique with
  respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations.

  Set a list of variables to some specified value for the new records

  - LHS refer to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, an expression or NA.

  For example:

        set_values_to = exprs(
          PARAMCD = "WOBS",
          PARAM = "Worst Observations"
        )

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   `NULL`

- keep_source_vars:

  Variables to keep from the source dataset

  For each event the specified variables are kept from the selected
  observations. The variables specified for `by_vars` and created by
  `set_values_to` are always kept. The `keep_source_vars` field of the
  event will take precedence over the value of the `keep_source_vars`
  argument.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `exprs(everything())`

## Value

The input dataset with the best or worst observation of each by group
added as new observations.

## Details

1.  For each event select the observations to consider:

    1.  If the event is of class `event`, the observations of the source
        dataset are restricted by `condition` and then the first or last
        (`mode`) observation per by group (`by_vars`) is selected.

        If the event is of class `event_joined`,
        [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
        is called to select the observations.

    2.  The variables specified by the `set_values_to` field of the
        event are added to the selected observations.

    3.  The variable specified for `tmp_event_nr_var` is added and set
        to the number of the event.

    4.  Only the variables specified for the `keep_source_vars` field of
        the event, and the by variables (`by_vars`) and the variables
        created by `set_values_to` are kept. If
        `keep_source_vars = NULL` is used for an event in
        `derive_extreme_event()` the value of the `keep_source_vars`
        argument of `derive_extreme_event()` is used.

2.  All selected observations are bound together.

3.  For each group (with respect to the variables specified for the
    `by_vars` parameter) the first or last observation (with respect to
    the order specified for the `order` parameter and the mode specified
    for the `mode` parameter) is selected.

4.  The variables specified by the `set_values_to` parameter are added
    to the selected observations.

5.  The observations are added to input dataset.

**Note:** This function creates temporary datasets which may be much
bigger than the input datasets. If this causes memory issues, please try
setting the admiral option `save_memory` to `TRUE` (see
[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)).
This reduces the memory consumption but increases the run-time.

## See also

[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)

## Examples

### Add a new record for the worst observation using [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md) objects

For each subject, the observation containing the worst sleeping problem
(if any exist) should be identified and added as a new record, retaining
all variables from the original observation. If multiple occurrences of
the worst sleeping problem occur, or no sleeping problems, then take the
observation occurring at the latest day.

- The groups for which new records are added are specified by the
  `by_vars` argument. Here for each *subject* a record should be added.
  Thus `by_vars = exprs(STUDYID, USUBJID)` is specified.

- The sets of possible sleeping problems are passed through the `events`
  argument as
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  objects. Each event contains a `condition` which may or may not be
  satisfied by each record (or possibly a group of records) within the
  input dataset `dataset`. Summary functions such as
  [`any()`](https://rdrr.io/r/base/any.html) and
  [`all()`](https://rdrr.io/r/base/all.html) are often handy to use
  within conditions, as is done here for the third event, which checks
  that the subject had no sleeping issues. The final event uses a
  catch-all `condition = TRUE` to ensure all subjects have a new record
  derived. Note that in this example, as no condition involves analysis
  of **cross-comparison values of within records**, it is sufficient to
  use
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  objects rather than
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md) -
  see the next example for a more complex condition.

- If any subject has one or more records satisfying the conditions from
  events, we can select just one record using the `order` argument. In
  this example, the first argument passed to `order` is `event_nr`,
  which is a temporary variable created through the `tmp_event_nr_var`
  argument, which numbers the events consecutively. Since
  `mode = "first"`, we only consider the first event for which a
  condition is satisfied. Within that event, we consider only the
  observation with the latest day, because the second argument for the
  order is `desc(ADY)`.

- Once a record is identified as satisfying an event's condition, a new
  observation is created by the following process:

  1.  the selected record is copied,

  2.  the variables specified in the event's `set_values_to` (here,
      `AVAL` and `AVALC`) are created/updated,

  3.  the variables specified in `keep_source_vars` (here, `ADY` does
      due to the use of the tidyselect expression
      [`everything()`](https://tidyselect.r-lib.org/reference/everything.html))
      (plus `by_vars` and the variables from `set_values_to`) are kept,

  4.  the variables specified in the global `set_values_to` (here,
      `PARAM` and `PARAMCD`) are created/updated.

    library(tibble, warn.conflicts = FALSE)
    library(dplyr, warn.conflicts = FALSE)
    library(lubridate, warn.conflicts = FALSE)

    adqs1 <- tribble(
      ~USUBJID, ~PARAMCD,         ~AVALC,        ~ADY,
      "1",      "NO SLEEP",       "N",              1,
      "1",      "WAKE UP 3X",     "N",              2,
      "2",      "NO SLEEP",       "N",              1,
      "2",      "WAKE UP 3X",     "Y",              2,
      "2",      "WAKE UP 3X",     "Y",              3,
      "3",      "NO SLEEP",       NA_character_,    1
    ) %>%
    mutate(STUDYID = "AB42")

    derive_extreme_event(
      adqs1,
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
          set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
        ),
        event(
          condition = PARAMCD == "WAKE UP 3X" & AVALC == "Y",
          set_values_to = exprs(AVALC = "Waking up three times", AVAL = 2)
        ),
        event(
          condition = all(AVALC == "N"),
          set_values_to = exprs(AVALC = "No sleeping problems", AVAL = 3)
        ),
        event(
          condition = TRUE,
          set_values_to = exprs(AVALC = "Missing", AVAL = 99)
        )
      ),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, desc(ADY)),
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "WSP",
        PARAM = "Worst Sleeping Problem"
      ),
      keep_source_vars = exprs(everything())
    ) %>%
    select(-STUDYID)
    #> # A tibble: 9 × 6
    #>   USUBJID PARAMCD    AVALC                   ADY  AVAL PARAM
    #>   <chr>   <chr>      <chr>                 <dbl> <dbl> <chr>
    #> 1 1       NO SLEEP   N                         1    NA <NA>
    #> 2 1       WAKE UP 3X N                         2    NA <NA>
    #> 3 2       NO SLEEP   N                         1    NA <NA>
    #> 4 2       WAKE UP 3X Y                         2    NA <NA>
    #> 5 2       WAKE UP 3X Y                         3    NA <NA>
    #> 6 3       NO SLEEP   <NA>                      1    NA <NA>
    #> 7 1       WSP        No sleeping problems      2     3 Worst Sleeping Problem
    #> 8 2       WSP        Waking up three times     3     2 Worst Sleeping Problem
    #> 9 3       WSP        Missing                   1    99 Worst Sleeping Problem

### Events based on comparison across records ([`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md))

We'll now extend the above example. Specifically, we consider a new
possible worst sleeping problem, namely if a subject experiences no
sleep on consecutive days.

- The "consecutive days" portion of the condition requires records to be
  compared with each other. This is done by using an
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  object, specifically by passing `dataset_name = adqs2` to it so that
  the `adqs2` dataset is joined onto itself. The `condition` now checks
  for two no sleep records, and crucially compares the `ADY` values to
  see if they differ by one day. The `.join` syntax distinguishes
  between the `ADY` value of the parent and joined datasets. As the
  condition involves `AVALC`, `PARAMCD` and `ADY`, we specify these
  variables with `join_vars`, and finally, because we wish to compare
  all records with each other, we select `join_type = "all"`.

    adqs2 <- tribble(
       ~USUBJID, ~PARAMCD,     ~AVALC, ~ADY,
       "4",      "WAKE UP",    "N",    1,
       "4",      "NO SLEEP",   "Y",    2,
       "4",      "NO SLEEP",   "Y",    3,
       "5",      "NO SLEEP",   "N",    1,
       "5",      "NO SLEEP",   "Y",    2,
       "5",      "WAKE UP 3X", "Y",    3,
       "5",      "NO SLEEP",   "Y",    4
    ) %>%
    mutate(STUDYID = "AB42")

    derive_extreme_event(
      adqs2,
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event_joined(
          join_vars = exprs(AVALC, PARAMCD, ADY),
          join_type = "all",
          condition = PARAMCD == "NO SLEEP" & AVALC == "Y" &
            PARAMCD.join == "NO SLEEP" & AVALC.join == "Y" &
            ADY == ADY.join + 1,
          set_values_to = exprs(AVALC = "No sleep two nights in a row", AVAL = 0)
        ),
        event(
          condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
          set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
        ),
        event(
          condition = PARAMCD == "WAKE UP 3X" & AVALC == "Y",
          set_values_to = exprs(AVALC = "Waking up three times", AVAL = 2)
        ),
        event(
          condition = all(AVALC == "N"),
          set_values_to = exprs(
            AVALC = "No sleeping problems", AVAL = 3
          )
        ),
        event(
          condition = TRUE,
          set_values_to = exprs(AVALC = "Missing", AVAL = 99)
        )
      ),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, desc(ADY)),
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "WSP",
        PARAM = "Worst Sleeping Problem"
      ),
      keep_source_vars = exprs(everything())
    ) %>%
    select(-STUDYID)
    #> # A tibble: 9 × 6
    #>   USUBJID PARAMCD    AVALC                          ADY  AVAL PARAM
    #>   <chr>   <chr>      <chr>                        <dbl> <dbl> <chr>
    #> 1 4       WAKE UP    N                                1    NA <NA>
    #> 2 4       NO SLEEP   Y                                2    NA <NA>
    #> 3 4       NO SLEEP   Y                                3    NA <NA>
    #> 4 5       NO SLEEP   N                                1    NA <NA>
    #> 5 5       NO SLEEP   Y                                2    NA <NA>
    #> 6 5       WAKE UP 3X Y                                3    NA <NA>
    #> 7 5       NO SLEEP   Y                                4    NA <NA>
    #> 8 4       WSP        No sleep two nights in a row     3     0 Worst Sleeping Pr…
    #> 9 5       WSP        No sleep                         4     1 Worst Sleeping Pr…

### Specifying different arguments across [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md) objects

Here we consider a Hy's Law use case. We are interested in knowing
whether a subject's Alkaline Phosphatase has ever been above twice the
upper limit of normal range. If so, i.e. if `CRIT1FL` is `Y`, we are
interested in the record for the first time this occurs, and if not, we
wish to retain the last record. As such, for this case now we need to
vary our usage of the `mode` argument dependent on the
[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).

- In first
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
  since we simply seek the first time that `CRIT1FL` is `"Y"`, it's
  enough to specify the `condition`, because we inherit `order` and
  `mode` from the main `derive_extreme_event()` call here which will
  automatically select the first occurrence by `AVISITN`.

- In the second
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
  we select the last record among the full set of records where
  `CRIT1FL` are all `"N"` by additionally specifying `mode = "last"`
  within the
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).

- Note now the usage of `keep_source_vars = exprs(AVISITN)` rather than
  [`everything()`](https://tidyselect.r-lib.org/reference/everything.html)
  as in the previous example. This is done to ensure `CRIT1` and
  `CRIT1FL` are not populated for the new records.

    adhy <- tribble(
      ~USUBJID, ~AVISITN,              ~CRIT1, ~CRIT1FL,
      "1",             1, "ALT > 2 times ULN", "N",
      "1",             2, "ALT > 2 times ULN", "N",
      "2",             1, "ALT > 2 times ULN", "N",
      "2",             2, "ALT > 2 times ULN", "Y",
      "2",             3, "ALT > 2 times ULN", "N",
      "2",             4, "ALT > 2 times ULN", "Y"
    ) %>%
      mutate(
        PARAMCD = "ALT",
        PARAM = "ALT (U/L)",
        STUDYID = "AB42"
      )

    derive_extreme_event(
      adhy,
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          condition = CRIT1FL == "Y",
          set_values_to = exprs(AVALC = "Y")
        ),
        event(
          condition = CRIT1FL == "N",
          mode = "last",
          set_values_to = exprs(AVALC = "N")
        )
      ),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, AVISITN),
      mode = "first",
      keep_source_vars = exprs(AVISITN),
      set_values_to = exprs(
        PARAMCD = "ALT2",
        PARAM = "ALT > 2 times ULN"
      )
    ) %>%
      select(-STUDYID)
    #> # A tibble: 8 × 7
    #>   USUBJID AVISITN CRIT1             CRIT1FL PARAMCD PARAM             AVALC
    #>   <chr>     <dbl> <chr>             <chr>   <chr>   <chr>             <chr>
    #> 1 1             1 ALT > 2 times ULN N       ALT     ALT (U/L)         <NA>
    #> 2 1             2 ALT > 2 times ULN N       ALT     ALT (U/L)         <NA>
    #> 3 2             1 ALT > 2 times ULN N       ALT     ALT (U/L)         <NA>
    #> 4 2             2 ALT > 2 times ULN Y       ALT     ALT (U/L)         <NA>
    #> 5 2             3 ALT > 2 times ULN N       ALT     ALT (U/L)         <NA>
    #> 6 2             4 ALT > 2 times ULN Y       ALT     ALT (U/L)         <NA>
    #> 7 1             2 <NA>              <NA>    ALT2    ALT > 2 times ULN N
    #> 8 2             2 <NA>              <NA>    ALT2    ALT > 2 times ULN Y    

### A more complex example: Confirmed Best Overall Response (`first/last_cond_upper`, `join_type`, `source_datasets`)

The final example showcases a use of `derive_extreme_event()` to
calculate the Confirmed Best Overall Response (CBOR) in an `ADRS`
dataset, as is common in many oncology trials. This example builds on
all the previous ones and thus assumes a baseline level of confidence
with `derive_extreme_event()`.

The following `ADSL` and `ADRS` datasets will be used throughout:

    adsl <- tribble(
      ~USUBJID, ~TRTSDTC,
      "1",      "2020-01-01",
      "2",      "2019-12-12",
      "3",      "2019-11-11",
      "4",      "2019-12-30",
      "5",      "2020-01-01",
      "6",      "2020-02-02",
      "7",      "2020-02-02",
      "8",      "2020-02-01"
    ) %>%
    mutate(
      TRTSDT = ymd(TRTSDTC),
      STUDYID = "AB42"
    )

    adrs <- tribble(
      ~USUBJID, ~ADTC,        ~AVALC,
      "1",      "2020-01-01", "PR",
      "1",      "2020-02-01", "CR",
      "1",      "2020-02-16", "NE",
      "1",      "2020-03-01", "CR",
      "1",      "2020-04-01", "SD",
      "2",      "2020-01-01", "SD",
      "2",      "2020-02-01", "PR",
      "2",      "2020-03-01", "SD",
      "2",      "2020-03-13", "CR",
      "4",      "2020-01-01", "PR",
      "4",      "2020-03-01", "NE",
      "4",      "2020-04-01", "NE",
      "4",      "2020-05-01", "PR",
      "5",      "2020-01-01", "PR",
      "5",      "2020-01-10", "PR",
      "5",      "2020-01-20", "PR",
      "6",      "2020-02-06", "PR",
      "6",      "2020-02-16", "CR",
      "6",      "2020-03-30", "PR",
      "7",      "2020-02-06", "PR",
      "7",      "2020-02-16", "CR",
      "7",      "2020-04-01", "NE",
      "8",      "2020-02-16", "PD"
    ) %>%
      mutate(
        ADT = ymd(ADTC),
        STUDYID = "AB42",
        PARAMCD = "OVR",
        PARAM = "Overall Response by Investigator"
      ) %>%
      derive_vars_merged(
        dataset_add = adsl,
        by_vars = exprs(STUDYID, USUBJID),
        new_vars = exprs(TRTSDT)
      )

Since the CBOR derivation contains multiple complex parts, it's
convenient to make use of the `description` argument within each event
object to describe what condition is being checked.

- For the Confirmed Response (CR), for each `"CR"` record in the
  original `ADRS` dataset that will be identified by the first part of
  the `condition` argument (`AVALC == "CR"`), we need to use the
  `first_cond_upper` argument to limit the group of observations to
  consider alongside it. Namely, we need to look up to and including the
  second CR (`AVALC.join == "CR"`) over 28 days from the first one
  (`ADT.join >= ADT + 28`). The observations satisfying
  `first_cond_upper` then form part of our "join group", meaning that
  the remaining portions of `condition` which reference joined variables
  are limited to this group. In particular, within `condition` we use
  [`all()`](https://rdrr.io/r/base/all.html) to check that all
  observations are either `"CR"` or `"NE"`, and
  [`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md)
  to ensure at most one is `"NE"`.

  Note that the selection of `join_type = "after"` is critical here, due
  to the fact that the restriction implied by `join_type` is applied
  before the one implied by `first_cond_upper`. Picking the first
  subject (who was correctly identified as a confirmed responder) as an
  example, selecting `join_type = "all"` instead of `"after"` would mean
  the first `"PR"` record from `"2020-01-01"` would also be considered
  when evaluating the `all(AVALC.join %in% c("CR", "NE"))` portion of
  `condition`. In turn, the condition would not be satisfied anymore,
  and in this case, following the later event logic shows the subject
  would be considered a partial responder instead.

- The Partial Response (PR), is very similar; with the difference being
  that the first portion of `condition` now references `"PR"` and
  `first_cond_upper` accepts a confirmatory `"PR"` or `"CR"` 28 days
  later. Note that now we must add `"PR"` as an option within the
  [`all()`](https://rdrr.io/r/base/all.html) condition to account for
  confirmatory `"PR"`s.

- The Stable Disease (SD), Progressive Disease (PD) and Not Evaluable
  (NE) events are simpler and just require
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  calls.

- Finally, we use a catch-all
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  with `condition = TRUE` and `dataset_name = "adsl"` to identify those
  subjects who do not appear in `ADRS` and list their CBOR as
  `"MISSING"`. Note here the fact that `dataset_name` is set to
  `"adsl"`, which is a new source dataset. As such it's important in the
  main `derive_extreme_event()` call to list `adsl` as another source
  dataset with `source_datasets = list(adsl = adsl)`.

    derive_extreme_event(
      adrs,
      by_vars = exprs(STUDYID, USUBJID),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, ADT),
      mode = "first",
      source_datasets = list(adsl = adsl),
      events = list(
        event_joined(
          description = paste(
            "CR needs to be confirmed by a second CR at least 28 days later",
            "at most one NE is acceptable between the two assessments"
          ),
          join_vars = exprs(AVALC, ADT),
          join_type = "after",
          first_cond_upper = AVALC.join == "CR" & ADT.join >= ADT + 28,
          condition = AVALC == "CR" &
            all(AVALC.join %in% c("CR", "NE")) &
            count_vals(var = AVALC.join, val = "NE") <= 1,
          set_values_to = exprs(AVALC = "CR")
        ),
        event_joined(
          description = paste(
            "PR needs to be confirmed by a second CR or PR at least 28 days later,",
            "at most one NE is acceptable between the two assessments"
          ),
          join_vars = exprs(AVALC, ADT),
          join_type = "after",
          first_cond_upper = AVALC.join %in% c("CR", "PR") & ADT.join >= ADT + 28,
          condition = AVALC == "PR" &
            all(AVALC.join %in% c("CR", "PR", "NE")) &
            count_vals(var = AVALC.join, val = "NE") <= 1,
          set_values_to = exprs(AVALC = "PR")
        ),
        event(
          description = paste(
            "CR, PR, or SD are considered as SD if occurring at least 28",
            "after treatment start"
          ),
          condition = AVALC %in% c("CR", "PR", "SD") & ADT >= TRTSDT + 28,
          set_values_to = exprs(AVALC = "SD")
        ),
        event(
          condition = AVALC == "PD",
          set_values_to = exprs(AVALC = "PD")
        ),
        event(
          condition = AVALC %in% c("CR", "PR", "SD", "NE"),
          set_values_to = exprs(AVALC = "NE")
        ),
        event(
          description = "Set response to MISSING for patients without records in ADRS",
          dataset_name = "adsl",
          condition = TRUE,
          set_values_to = exprs(AVALC = "MISSING"),
          keep_source_vars = exprs(TRTSDT)
        )
      ),
      set_values_to = exprs(
        PARAMCD = "CBOR",
        PARAM = "Best Confirmed Overall Response by Investigator"
      )
    ) %>%
      filter(PARAMCD == "CBOR") %>%
      select(-STUDYID, -ADTC)
    #> # A tibble: 8 × 6
    #>   USUBJID AVALC   ADT        PARAMCD PARAM                            TRTSDT
    #>   <chr>   <chr>   <date>     <chr>   <chr>                            <date>
    #> 1 1       CR      2020-02-01 CBOR    Best Confirmed Overall Response… 2020-01-01
    #> 2 2       SD      2020-02-01 CBOR    Best Confirmed Overall Response… 2019-12-12
    #> 3 3       MISSING NA         CBOR    Best Confirmed Overall Response… 2019-11-11
    #> 4 4       SD      2020-05-01 CBOR    Best Confirmed Overall Response… 2019-12-30
    #> 5 5       NE      2020-01-01 CBOR    Best Confirmed Overall Response… 2020-01-01
    #> 6 6       PR      2020-02-06 CBOR    Best Confirmed Overall Response… 2020-02-02
    #> 7 7       NE      2020-02-06 CBOR    Best Confirmed Overall Response… 2020-02-02
    #> 8 8       PD      2020-02-16 CBOR    Best Confirmed Overall Response… 2020-02-01

### Further examples

Equivalent examples for using the`check_type` argument can be found in
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md).
