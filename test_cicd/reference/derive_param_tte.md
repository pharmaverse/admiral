# Derive a Time-to-Event Parameter

Add a time-to-event parameter to the input dataset.

## Usage

``` r
derive_param_tte(
  dataset = NULL,
  dataset_adsl,
  source_datasets,
  by_vars = NULL,
  start_date = TRTSDT,
  event_conditions,
  censor_conditions,
  create_datetime = FALSE,
  set_values_to,
  subject_keys = get_admiral_option("subject_keys"),
  check_type = "warning"
)
```

## Arguments

- dataset:

  Input dataset

  `PARAMCD` is expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- dataset_adsl:

  ADSL input dataset

  The variables specified for `start_date`, and `subject_keys` are
  expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- source_datasets:

  Source datasets

  A named list of datasets is expected. The `dataset_name` field of
  [`tte_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/tte_source.md)
  refers to the dataset provided in the list.

  Permitted values

  :   named list of datasets, e.g., `list(adsl = adsl, ae = ae)`

  Default value

  :   none

- by_vars:

  By variables

  If the parameter is specified, separate time to event parameters are
  derived for each by group.

  The by variables must be in at least one of the source datasets. Each
  source dataset must contain either all by variables or none of the by
  variables.

  The by variables are not included in the output dataset.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- start_date:

  Time to event origin date

  The variable `STARTDT` is set to the specified date. The value is
  taken from the ADSL dataset.

  If the event or censoring date is before the origin date, `ADT` is set
  to the origin date.

  Permitted values

  :   a date or datetime variable

  Default value

  :   `TRTSDT`

- event_conditions:

  Sources and conditions defining events

  A list of
  [`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)
  objects is expected.

  Permitted values

  :   a list of source objects, e.g., `list(pd, death)`

  Default value

  :   none

- censor_conditions:

  Sources and conditions defining censorings

  A list of
  [`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
  objects is expected.

  Permitted values

  :   a list of source objects, e.g., `list(pd, death)`

  Default value

  :   none

- create_datetime:

  Create datetime variables?

  If set to `TRUE`, variables `ADTM` and `STARTDTM` are created.
  Otherwise, variables `ADT` and `STARTDT` are created.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- set_values_to:

  Variables to set

  A named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining the variables to be set for the new parameter, e.g.
  `exprs(PARAMCD = "OS", PARAM = "Overall Survival")` is expected. The
  values must be symbols, character strings, numeric values,
  expressions, or `NA`.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of symbols created using
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  is expected.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `get_admiral_option("subject_keys")`

- check_type:

  Check uniqueness

  If `"warning"`, `"message"`, or `"error"` is specified, the specified
  message is issued if the observations of the source datasets are not
  unique with respect to the by variables and the date and order
  specified in the
  [`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)
  and
  [`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
  objects.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

The input dataset with the new parameter added

## Details

The following steps are performed to create the observations of the new
parameter:

**Deriving the events:**

1.  For each event source dataset the observations as specified by the
    `filter` element are selected. Then for each subject the first
    observation (with respect to `date` and `order`) is selected.

2.  The `ADT` variable is set to the variable specified by the `date`
    element. If the date variable is a datetime variable, only the
    datepart is copied.

3.  The `CNSR` variable is added and set to the `censor` element.

4.  The variables specified by the `set_values_to` element are added.

5.  The selected observations of all event source datasets are combined
    into a single dataset.

6.  For each subject the first observation (with respect to the
    `ADT`/`ADTM` variable) from the single dataset is selected. If there
    is more than one event with the same date, the first event with
    respect to the order of events in `event_conditions` is selected.

**Deriving the censoring observations:**

1.  For each censoring source dataset the observations as specified by
    the `filter` element are selected. Then for each subject the last
    observation (with respect to `date` and `order`) is selected.

2.  The `ADT` variable is set to the variable specified by the `date`
    element. If the date variable is a datetime variable, only the
    datepart is copied.

3.  The `CNSR` variable is added and set to the `censor` element.

4.  The variables specified by the `set_values_to` element are added.

5.  The selected observations of all censoring source datasets are
    combined into a single dataset.

6.  For each subject the last observation (with respect to the
    `ADT`/`ADTM` variable) from the single dataset is selected. If there
    is more than one censoring with the same date, the last censoring
    with respect to the order of censorings in `censor_conditions` is
    selected.

For each subject (as defined by the `subject_keys` parameter) an
observation is selected. If an event is available, the event observation
is selected. Otherwise the censoring observation is selected.

Finally:

1.  The variable specified for `start_date` is joined from the ADSL
    dataset. Only subjects in both datasets are kept, i.e., subjects
    with both an event or censoring and an observation in
    `dataset_adsl`.

2.  The variables as defined by the `set_values_to` parameter are added.

3.  The `ADT`/`ADTM` variable is set to the maximum of `ADT`/`ADTM` and
    `STARTDT`/`STARTDTM` (depending on the `create_datetime` parameter).

4.  The new observations are added to the output dataset.

## See also

[`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)

## Examples

### Add a basic time to event parameter

For each subject the time to first adverse event should be created as a
parameter.

- The event source object is created using
  [`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)
  and the date is set to adverse event start date.

- The censor source object is created using
  [`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
  and the date is set to end of study date.

- The event and censor source objects are then passed to
  `derive_param_tte()` to derive the time to event parameter with the
  provided parameter descriptions (`PARAMCD` and `PARAM`).

- Note the values of the censor variable (`CNSR`) that are derived
  below, where the first subject has an event and the second does not.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)
    library(lubridate, warn.conflicts = FALSE)

    adsl <- tribble(
      ~USUBJID, ~TRTSDT,           ~EOSDT,            ~NEWDRGDT,
      "01",     ymd("2020-12-06"), ymd("2021-03-06"), NA,
      "02",     ymd("2021-01-16"), ymd("2021-02-03"), ymd("2021-01-03")
    ) %>%
      mutate(STUDYID = "AB42")

    adae <- tribble(
      ~USUBJID, ~ASTDT,            ~AESEQ, ~AEDECOD,
      "01",     ymd("2021-01-03"),      1, "Flu",
      "01",     ymd("2021-03-04"),      2, "Cough",
      "01",     ymd("2021-03-05"),      3, "Cough"
    ) %>%
      mutate(STUDYID = "AB42")

    ttae <- event_source(
      dataset_name = "adae",
      date = ASTDT,
      set_values_to = exprs(
        EVNTDESC = "AE",
        SRCDOM = "ADAE",
        SRCVAR = "ASTDT",
        SRCSEQ = AESEQ
      )
    )

    eos <- censor_source(
      dataset_name = "adsl",
      date = EOSDT,
      set_values_to = exprs(
        EVNTDESC = "END OF STUDY",
        SRCDOM = "ADSL",
        SRCVAR = "EOSDT"
      )
    )

    derive_param_tte(
      dataset_adsl = adsl,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, adae = adae),
      set_values_to = exprs(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 2 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAE    Time to First Adverse Event 2021-01-03     0      1
    #> 2 02      2021-01-16 TTAE    Time to First Adverse Event 2021-02-03     1     NA

### Adding a by variable (`by_vars`)

By variables can be added using the `by_vars` argument, e.g., now for
each subject the time to first occurrence of each adverse event
preferred term (`AEDECOD`) should be created as parameters.

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, adae = adae),
      set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event")
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 4 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAE1   Time to First Cough Advers… 2021-03-04     0      2
    #> 2 01      2020-12-06 TTAE2   Time to First Flu Adverse … 2021-01-03     0      1
    #> 3 02      2021-01-16 TTAE1   Time to First Cough Advers… 2021-02-03     1     NA
    #> 4 02      2021-01-16 TTAE2   Time to First Flu Adverse … 2021-02-03     1     NA

### Handling duplicates (`check_type`)

The source records are checked regarding duplicates with respect to the
by variables and the date and order specified in the source objects. By
default, a warning is issued if any duplicates are found. Note here how
after creating a new adverse event dataset containing a duplicate date
for `"Cough"`, it was then passed to the function using the
`source_datasets` argument - where you see below `adae = adae_dup`.

    adae_dup <- tribble(
      ~USUBJID, ~ASTDT,            ~AESEQ, ~AEDECOD, ~AESER,
      "01",     ymd("2021-01-03"),      1, "Flu",    "Y",
      "01",     ymd("2021-03-04"),      2, "Cough",  "N",
      "01",     ymd("2021-03-04"),      3, "Cough",  "Y"
    ) %>%
      mutate(STUDYID = "AB42")

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      start_date = TRTSDT,
      source_datasets = list(adsl = adsl, adae = adae_dup),
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event")
      )
    )
    #> # A tibble: 4 × 11
    #>   USUBJID STUDYID EVNTDESC     SRCDOM SRCVAR SRCSEQ  CNSR ADT        STARTDT
    #>   <chr>   <chr>   <chr>        <chr>  <chr>   <dbl> <int> <date>     <date>
    #> 1 01      AB42    AE           ADAE   ASTDT       2     0 2021-03-04 2020-12-06
    #> 2 01      AB42    AE           ADAE   ASTDT       1     0 2021-01-03 2020-12-06
    #> 3 02      AB42    END OF STUDY ADSL   EOSDT      NA     1 2021-02-03 2021-01-16
    #> 4 02      AB42    END OF STUDY ADSL   EOSDT      NA     1 2021-02-03 2021-01-16
    #> # i 2 more variables: PARAMCD <chr>, PARAM <chr>
    #> Warning: Dataset "adae" contains duplicate records with respect to `STUDYID`, `USUBJID`,
    #> `AEDECOD`, and `ASTDT`
    #> i Run `admiral::get_duplicates_dataset()` to access the duplicate records

For investigating the issue, the dataset of the duplicate source records
can be obtained by calling
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_duplicates_dataset.md):

    get_duplicates_dataset()
    #> Duplicate records with respect to `STUDYID`, `USUBJID`, `AEDECOD`, and `ASTDT`.
    #> # A tibble: 2 × 6
    #>   STUDYID USUBJID AEDECOD ASTDT      AESEQ AESER
    #> * <chr>   <chr>   <chr>   <date>     <dbl> <chr>
    #> 1 AB42    01      Cough   2021-03-04     2 N
    #> 2 AB42    01      Cough   2021-03-04     3 Y    

Common options to solve the issue:

- Restricting the source records by specifying/updating the `filter`
  argument in the
  [`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)/[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
  calls.

- Specifying additional variables for `order` in the
  [`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)/[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
  calls.

- Setting `check_type = "none"` in the `derive_param_tte()` call to
  ignore any duplicates.

In this example it does not have significant impact which record is
chosen as the dates are the same so the time to event derivation will be
the same, but it does impact `SRCSEQ` in the output dataset, so here the
second option is used. Note here how you can also define source objects
from within the `derive_param_tte()` function call itself.

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      start_date = TRTSDT,
      source_datasets = list(adsl = adsl, adae = adae_dup),
      event_conditions = list(event_source(
        dataset_name = "adae",
        date = ASTDT,
        set_values_to = exprs(
          EVNTDESC = "AE",
          SRCDOM = "ADAE",
          SRCVAR = "ASTDT",
          SRCSEQ = AESEQ
        ),
        order = exprs(AESEQ)
      )),
      censor_conditions = list(eos),
      set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event")
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 4 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAE1   Time to First Cough Advers… 2021-03-04     0      2
    #> 2 01      2020-12-06 TTAE2   Time to First Flu Adverse … 2021-01-03     0      1
    #> 3 02      2021-01-16 TTAE1   Time to First Cough Advers… 2021-02-03     1     NA
    #> 4 02      2021-01-16 TTAE2   Time to First Flu Adverse … 2021-02-03     1     NA

### Filtering source records (`filter`)

The first option from above could have been achieved using `filter`, for
example here only using serious adverse events.

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      start_date = TRTSDT,
      source_datasets = list(adsl = adsl, adae = adae_dup),
      event_conditions = list(event_source(
        dataset_name = "adae",
        filter = AESER == "Y",
        date = ASTDT,
        set_values_to = exprs(
          EVNTDESC = "Serious AE",
          SRCDOM = "ADAE",
          SRCVAR = "ASTDT",
          SRCSEQ = AESEQ
        )
      )),
      censor_conditions = list(eos),
      set_values_to = exprs(
        PARAMCD = paste0("TTSAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First Serious", AEDECOD, "Adverse Event")
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 4 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTSAE1  Time to First Serious Coug… 2021-03-04     0      3
    #> 2 01      2020-12-06 TTSAE2  Time to First Serious Flu … 2021-01-03     0      1
    #> 3 02      2021-01-16 TTSAE1  Time to First Serious Coug… 2021-02-03     1     NA
    #> 4 02      2021-01-16 TTSAE2  Time to First Serious Flu … 2021-02-03     1     NA

### Using multiple event/censor conditions (`event_conditions` /`censor_conditions`)

In the above examples, we only have a single event and single censor
condition. Here, we now consider multiple conditions for each passed
using `event_conditions` and `censor_conditions`.

For the event we are going to use first AE and additionally check a lab
condition, and for the censor we'll add in treatment start date in case
end of study date was ever missing.

    adlb <- tribble(
      ~USUBJID, ~ADT,              ~PARAMCD, ~ANRIND,
      "01",     ymd("2020-12-22"), "HGB",    "LOW"
    ) %>%
      mutate(STUDYID = "AB42")

    low_hgb <- event_source(
      dataset_name = "adlb",
      filter = PARAMCD == "HGB" & ANRIND == "LOW",
      date = ADT,
      set_values_to = exprs(
        EVNTDESC = "POSSIBLE ANEMIA",
        SRCDOM = "ADLB",
        SRCVAR = "ADT"
      )
    )

    trt_start <- censor_source(
      dataset_name = "adsl",
      date = TRTSDT,
      set_values_to = exprs(
        EVNTDESC = "TREATMENT START",
        SRCDOM = "ADSL",
        SRCVAR = "TRTSDT"
      )
    )

    derive_param_tte(
      dataset_adsl = adsl,
      event_conditions = list(ttae, low_hgb),
      censor_conditions = list(eos, trt_start),
      source_datasets = list(adsl = adsl, adae = adae, adlb = adlb),
      set_values_to = exprs(
        PARAMCD = "TTAELB",
        PARAM = "Time to First Adverse Event or Possible Anemia (Labs)"
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 2 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAELB  Time to First Adverse Even… 2020-12-22     0     NA
    #> 2 02      2021-01-16 TTAELB  Time to First Adverse Even… 2021-02-03     1     NA

Note above how the earliest event date is always taken and the latest
censor date.

### Using different censor values (`censor`) and censoring at earliest occurring censor condition

Within
[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md)
the value used to denote a censor can be changed from the default of
`1`.

In this example an extra censor is used for new drug date with the value
of `2`.

    newdrug <- censor_source(
      dataset_name = "adsl",
      date = NEWDRGDT,
      censor = 2,
      set_values_to = exprs(
        EVNTDESC = "NEW DRUG RECEIVED",
        SRCDOM = "ADSL",
        SRCVAR = "NEWDRGDT"
      )
    )

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      event_conditions = list(ttae),
      censor_conditions = list(eos, newdrug),
      source_datasets = list(adsl = adsl, adae = adae),
      set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event")
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 4 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAE1   Time to First Cough Advers… 2021-03-04     0      2
    #> 2 01      2020-12-06 TTAE2   Time to First Flu Adverse … 2021-01-03     0      1
    #> 3 02      2021-01-16 TTAE1   Time to First Cough Advers… 2021-02-03     1     NA
    #> 4 02      2021-01-16 TTAE2   Time to First Flu Adverse … 2021-02-03     1     NA

In this case the results are still the same, because as explained in the
above example the latest censor condition is always taken for those
without an event. For the second subject this is still the end of study
date.

So, if we wanted to instead censor here at the new drug date if subject
has one, then we would need to again use the `filter` argument, but this
time for a new end of study censor source object.

    eos_nonewdrug <- censor_source(
      dataset_name = "adsl",
      filter = is.na(NEWDRGDT),
      date = EOSDT,
      set_values_to = exprs(
        EVNTDESC = "END OF STUDY",
        SRCDOM = "ADSL",
        SRCVAR = "EOSDT"
      )
    )

    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = exprs(AEDECOD),
      event_conditions = list(ttae),
      censor_conditions = list(eos_nonewdrug, newdrug),
      source_datasets = list(adsl = adsl, adae = adae),
      set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event")
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 4 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                       ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                       <date>     <int>  <dbl>
    #> 1 01      2020-12-06 TTAE1   Time to First Cough Advers… 2021-03-04     0      2
    #> 2 01      2020-12-06 TTAE2   Time to First Flu Adverse … 2021-01-03     0      1
    #> 3 02      2021-01-16 TTAE1   Time to First Cough Advers… 2021-01-16     2     NA
    #> 4 02      2021-01-16 TTAE2   Time to First Flu Adverse … 2021-01-16     2     NA

### Overall survival time to event parameter

In oncology trials, this is commonly derived as time from randomization
date to death. For those without event, they are censored at the last
date they are known to be alive.

- The start date is set using `start_date` argument, now that we need to
  use different to the default.

- In this example, datetime was needed, which can be achieved by setting
  `create_datetime` argument to `TRUE`.

    adsl <- tribble(
      ~USUBJID, ~RANDDTM,                       ~LSALVDTM,                      ~DTHDTM,                        ~DTHFL,
      "01",     ymd_hms("2020-10-03 00:00:00"), ymd_hms("2022-12-15 23:59:59"), NA,                             NA,
      "02",     ymd_hms("2021-01-23 00:00:00"), ymd_hms("2021-02-03 19:45:59"), ymd_hms("2021-02-03 19:45:59"), "Y"
    ) %>%
      mutate(STUDYID = "AB42")

    # derive overall survival parameter
    death <- event_source(
      dataset_name = "adsl",
      filter = DTHFL == "Y",
      date = DTHDTM,
      set_values_to = exprs(
        EVNTDESC = "DEATH",
        SRCDOM = "ADSL",
        SRCVAR = "DTHDTM"
      )
    )

    last_alive <- censor_source(
      dataset_name = "adsl",
      date = LSALVDTM,
      set_values_to = exprs(
        EVNTDESC = "LAST DATE KNOWN ALIVE",
        SRCDOM = "ADSL",
        SRCVAR = "LSALVDTM"
      )
    )

    derive_param_tte(
      dataset_adsl = adsl,
      start_date = RANDDTM,
      event_conditions = list(death),
      censor_conditions = list(last_alive),
      create_datetime = TRUE,
      source_datasets = list(adsl = adsl),
      set_values_to = exprs(
        PARAMCD = "OS",
        PARAM = "Overall Survival"
      )
    ) %>%
      select(USUBJID, STARTDTM, PARAMCD, PARAM, ADTM, CNSR)
    #> # A tibble: 2 × 6
    #>   USUBJID STARTDTM            PARAMCD PARAM            ADTM                 CNSR
    #>   <chr>   <dttm>              <chr>   <chr>            <dttm>              <int>
    #> 1 01      2020-10-03 00:00:00 OS      Overall Survival 2022-12-15 23:59:59     1
    #> 2 02      2021-01-23 00:00:00 OS      Overall Survival 2021-02-03 19:45:59     0

### Duration of response time to event parameter

In oncology trials, this is commonly derived as time from response until
progression or death, or if neither have occurred then censor at last
tumor assessment visit date. It is only relevant for subjects with a
response. Note how only observations for subjects in `dataset_adsl` have
the new parameter created, so see below how this is filtered only on
responders.

    adsl_resp <- tribble(
      ~USUBJID, ~DTHFL, ~DTHDT,            ~RSPDT,
      "01",     "Y",    ymd("2021-06-12"), ymd("2021-03-04"),
      "02",     "N",    NA,                NA,
      "03",     "Y",    ymd("2021-08-21"), NA,
      "04",     "N",    NA,                ymd("2021-04-14")
    ) %>%
      mutate(STUDYID = "AB42")

    adrs <- tribble(
      ~USUBJID, ~AVALC, ~ADT,              ~ASEQ,
      "01",     "SD",   ymd("2021-01-03"), 1,
      "01",     "PR",   ymd("2021-03-04"), 2,
      "01",     "PD",   ymd("2021-05-05"), 3,
      "02",     "PD",   ymd("2021-02-03"), 1,
      "04",     "SD",   ymd("2021-02-13"), 1,
      "04",     "PR",   ymd("2021-04-14"), 2,
      "04",     "CR",   ymd("2021-05-15"), 3
    ) %>%
      mutate(STUDYID = "AB42", PARAMCD = "OVR")

    pd <- event_source(
      dataset_name = "adrs",
      filter = AVALC == "PD",
      date = ADT,
      set_values_to = exprs(
        EVENTDESC = "PD",
        SRCDOM = "ADRS",
        SRCVAR = "ADTM",
        SRCSEQ = ASEQ
      )
    )

    death <- event_source(
      dataset_name = "adsl",
      filter = DTHFL == "Y",
      date = DTHDT,
      set_values_to = exprs(
        EVENTDESC = "DEATH",
        SRCDOM = "ADSL",
        SRCVAR = "DTHDT"
      )
    )

    last_visit <- censor_source(
      dataset_name = "adrs",
      date = ADT,
      set_values_to = exprs(
        EVENTDESC = "LAST TUMOR ASSESSMENT",
        SRCDOM = "ADRS",
        SRCVAR = "ADTM",
        SRCSEQ = ASEQ
      )
    )

    derive_param_tte(
      dataset_adsl = filter(adsl_resp, !is.na(RSPDT)),
      start_date = RSPDT,
      event_conditions = list(pd, death),
      censor_conditions = list(last_visit),
      source_datasets = list(adsl = adsl_resp, adrs = adrs),
      set_values_to = exprs(
        PARAMCD = "DURRSP",
        PARAM = "Duration of Response"
      )
    ) %>%
      select(USUBJID, STARTDT, PARAMCD, PARAM, ADT, CNSR, SRCSEQ)
    #> # A tibble: 2 × 7
    #>   USUBJID STARTDT    PARAMCD PARAM                ADT         CNSR SRCSEQ
    #>   <chr>   <date>     <chr>   <chr>                <date>     <int>  <dbl>
    #> 1 01      2021-03-04 DURRSP  Duration of Response 2021-05-05     0      3
    #> 2 04      2021-04-14 DURRSP  Duration of Response 2021-05-15     1      3

### Further examples

Further example usages of this function can be found in the
[`vignette("bds_tte")`](https:/pharmaverse.github.io/admiral/test_cicd/articles/bds_tte.md).
