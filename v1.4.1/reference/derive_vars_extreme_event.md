# Add the Worst or Best Observation for Each By Group as New Variables

Add the first available record from `events` for each by group as new
variables, all variables of the selected observation are kept. It can be
used for selecting the extreme observation from a series of user-defined
events.

## Usage

``` r
derive_vars_extreme_event(
  dataset,
  by_vars,
  events,
  tmp_event_nr_var = NULL,
  order,
  mode,
  source_datasets = NULL,
  check_type = "warning",
  new_vars
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- events:

  Conditions and new values defining events

  A list of
  [`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md)
  or
  [`event_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_joined.md)
  objects is expected. Only observations listed in the `events` are
  considered for deriving extreme event. If multiple records meet the
  filter `condition`, take the first record sorted by `order`. The data
  is grouped by `by_vars`, i.e., summary functions like
  [`all()`](https://rdrr.io/r/base/all.html) or
  [`any()`](https://rdrr.io/r/base/any.html) can be used in `condition`.

  For
  [`event_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_joined.md)
  events the observations are selected by calling
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_joined.md).
  The `condition` field is passed to the `filter_join` argument.

  Default value

  :   none

- tmp_event_nr_var:

  Temporary event number variable

  The specified variable is added to all source datasets and is set to
  the number of the event before selecting the records of the event.

  It can be used in `order` to determine which record should be used if
  records from more than one event are selected.

  The variable is not included in the output dataset.

  Default value

  :   `NULL`

- order:

  Sort order

  If a particular event from `events` has more than one observation,
  within the event and by group, the records are ordered by the
  specified order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))`

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
  [`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md)
  and
  [`event_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_joined.md)
  refers to the dataset provided in the list.

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

- new_vars:

  Variables to add

  The specified variables from the events are added to the output
  dataset. Variables can be renamed by naming the element, i.e.,
  `new_vars = exprs(<new name> = <old name>)`.

  Default value

  :   none

## Value

The input dataset with the best or worst observation of each by group
added as new variables.

## Details

1.  For each event select the observations to consider:

    1.  If the event is of class `event`, the observations of the source
        dataset are restricted by `condition` and then the first or last
        (`mode`) observation per by group (`by_vars`) is selected.

        If the event is of class `event_joined`,
        [`filter_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_joined.md)
        is called to select the observations.

    2.  The variables specified by the `set_values_to` field of the
        event are added to the selected observations.

    3.  The variable specified for `tmp_event_nr_var` is added and set
        to the number of the event.

2.  All selected observations are bound together.

3.  For each group (with respect to the variables specified for the
    `by_vars` parameter) the first or last observation (with respect to
    the order specified for the `order` parameter and the mode specified
    for the `mode` parameter) is selected.

4.  The variables specified by the `new_vars` parameter are added to the
    selected observations.

5.  The variables are added to input dataset.

## See also

[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_joined.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_event.md)

ADSL Functions that returns variable appended to dataset:
[`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_age_years.md),
[`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_aage.md),
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md)

## Examples

``` r
library(tibble)
library(dplyr)
library(lubridate)

adsl <- tribble(
  ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
  "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
  "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""),
  "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""),
  "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""),
  "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
)

lb <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
  "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
  "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
  "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
  "PILOT01",    "LB", "01-1133",    304, "2013-05-01T10:13",
  "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
  "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
  "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
  "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
  "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
  "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
) %>%
  mutate(
    ADT = convert_dtc_to_dt(LBDTC)
  )

derive_vars_extreme_event(
  adsl,
  by_vars = exprs(STUDYID, USUBJID),
  events = list(
    event(
      dataset_name = "adsl",
      condition = !is.na(DTHDT),
      set_values_to = exprs(LSTALVDT = DTHDT, DTHFL = "Y")
    ),
    event(
      dataset_name = "lb",
      condition = !is.na(ADT),
      order = exprs(ADT),
      mode = "last",
      set_values_to = exprs(LSTALVDT = ADT, DTHFL = "N")
    ),
    event(
      dataset_name = "adsl",
      condition = !is.na(TRTEDT),
      order = exprs(TRTEDT),
      mode = "last",
      set_values_to = exprs(LSTALVDT = TRTEDT, DTHFL = "N")
    )
  ),
  source_datasets = list(adsl = adsl, lb = lb),
  tmp_event_nr_var = event_nr,
  order = exprs(LSTALVDT, event_nr),
  mode = "last",
  new_vars = exprs(LSTALVDT, DTHFL)
)
#> # A tibble: 5 × 6
#>   STUDYID USUBJID TRTEDT     DTHDT      LSTALVDT   DTHFL
#>   <chr>   <chr>   <date>     <date>     <date>     <chr>
#> 1 PILOT01 01-1130 2014-08-16 2014-09-13 2014-09-13 Y    
#> 2 PILOT01 01-1133 2013-04-28 NA         2013-05-01 N    
#> 3 PILOT01 01-1211 2013-01-12 NA         2013-01-12 N    
#> 4 PILOT01 09-1081 2014-04-27 NA         2014-05-10 N    
#> 5 PILOT01 09-1088 2014-10-09 2014-11-01 2014-11-01 Y    

# Derive DTHCAUS from AE and DS domain data
adsl <- tribble(
  ~STUDYID,  ~USUBJID,
  "STUDY01", "PAT01",
  "STUDY01", "PAT02",
  "STUDY01", "PAT03"
)
ae <- tribble(
  ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
  "STUDY01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04",
  "STUDY01", "PAT01", 13, "CARDIAC ARREST", "FATAL", "2021-04-03",
)

ds <- tribble(
  ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
  "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
  "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
  "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
  "STUDY01", "PAT03", 1, "DEATH", "POST STUDY REPORTING OF DEATH", "2022-03-03"
)

derive_vars_extreme_event(
  adsl,
  by_vars = exprs(STUDYID, USUBJID),
  events = list(
    event(
      dataset_name = "ae",
      condition = AEOUT == "FATAL",
      set_values_to = exprs(DTHCAUS = AEDECOD, DTHDT = convert_dtc_to_dt(AEDTHDTC)),
      order = exprs(DTHDT)
    ),
    event(
      dataset_name = "ds",
      condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
      set_values_to = exprs(DTHCAUS = DSTERM, DTHDT = convert_dtc_to_dt(DSSTDTC)),
      order = exprs(DTHDT)
    )
  ),
  source_datasets = list(ae = ae, ds = ds),
  tmp_event_nr_var = event_nr,
  order = exprs(DTHDT, event_nr),
  mode = "first",
  new_vars = exprs(DTHCAUS, DTHDT)
)
#> # A tibble: 3 × 4
#>   STUDYID USUBJID DTHCAUS                             DTHDT     
#>   <chr>   <chr>   <chr>                               <date>    
#> 1 STUDY01 PAT01   CARDIAC ARREST                      2021-04-03
#> 2 STUDY01 PAT02   DEATH DUE TO PROGRESSION OF DISEASE 2022-02-01
#> 3 STUDY01 PAT03   NA                                  NA        
```
