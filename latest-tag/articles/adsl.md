# Creating ADSL

## Introduction

This article describes creating an `ADSL` ADaM. Examples are currently
presented and tested using `DM`, `EX` , `AE`, `LB` and `DS` SDTM
domains. However, other domains could be used.

**Note:** *All examples assume CDISC SDTM and/or ADaM format as input
unless otherwise specified.*

## Programming Flow

- [Read in Data](#readdata)
- [Derive Period, Subperiod, and Phase Variables (e.g. `APxxSDT`,
  `APxxEDT`, …)](#periodvars)
- [Derive Treatment Variables (`TRT0xP`, `TRT0xA`)](#treatmentvar)
- [Derive/Impute Numeric Treatment Date/Time and Duration (`TRTSDT`,
  `TRTEDT`, `TRTDURD`)](#trtdatetime)
- [Derive Disposition Variables](#disposition)
  - [Disposition Dates (e.g. `EOSDT`)](#disposition_date)
  - [Disposition Status (e.g. `EOSTT`)](#disposition_status)
  - [Disposition Reason(s) (e.g. `DCSREAS`,
    `DCSREASP`)](#disposition_reason)
  - [Randomization Date (`RANDDT`)](#randomization_date)
- [Derive Birth Date and Analysis Age (`BRTHDT`, `AAGE`,
  `AAGEU`)](#aage)
- [Derive Death Variables](#death)
  - [Death Date (`DTHDT`)](#death_date)
  - [Cause of Death (`DTHCAUS`)](#death_cause)
  - [Duration Relative to Death](#death_other)
- [Derive Last Known Date Alive (`LSTALVDT`)](#lstalvdt)
- [Derive Groupings and Populations](#groupings)
  - [Grouping (e.g. `AGEGR1` or `REGION1`)](#groupings_ex)
  - [Population Flags (e.g. `SAFFL`)](#popflag)
- [Derive Other Variables](#other)
- [Add Labels and Attributes](#attributes)

### Read in Data

To start, all data frames needed for the creation of `ADSL` should be
read into the environment. This will be a company specific process. Some
of the data frames needed may be `DM`, `EX`, `DS`, `AE`, and `LB`.

For example purpose, the CDISC Pilot SDTM datasets—which are included in
[pharmaversesdtm](https://pharmaverse.github.io/pharmaversesdtm/)—are
used.

``` r
library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)

dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)
```

The `DM` domain is used as the basis for `ADSL`:

``` r
adsl <- dm %>%
  select(-DOMAIN)
```

### Derive Period, Subperiod, and Phase Variables (e.g. `APxxSDT`, `APxxEDT`, …)

See the [“Visit and Period Variables”
vignette](https:/pharmaverse.github.io/admiral/v1.4.1/articles/visits_periods.html#periods_adsl)
for more information.

If the variables are not derived based on a period reference dataset,
they may be derived at a later point of the flow. For example, phases
like “Treatment Phase” and “Follow up” could be derived based on
treatment start and end date.

### Derive Treatment Variables (`TRT0xP`, `TRT0xA`)

The mapping of the treatment variables is left to the ADaM programmer.
An example mapping for a study without periods may be:

``` r
adsl <- dm %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM)
```

For studies with periods see the [“Visit and Period Variables”
vignette](https:/pharmaverse.github.io/admiral/v1.4.1/articles/visits_periods.html#treatment_adsl).

### Derive/Impute Numeric Treatment Date/Time and Duration (`TRTSDTM`, `TRTEDTM`, `TRTDURD`)

The function
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
can be used to derive the treatment start and end date/times using the
`ex` domain. A pre-processing step for `ex` is required to convert the
variable `EXSTDTC` and `EXSTDTC` to datetime variables and impute
missing date or time components. Conversion and imputation is done by
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm.md).

Example calls:

``` r
# Impute start and end time of exposure to first and last respectively,
# Do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )
```

This call returns the original data frame with the column `TRTSDTM`,
`TRTSTMF`, `TRTEDTM`, and `TRTETMF` added. Exposure observations with
incomplete date and zero doses of non placebo treatments are ignored.
Missing time parts are imputed as first or last for start and end date
respectively.

The datetime variables returned can be converted to dates using the
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm_to_dt.md)
function.

``` r
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))
```

Now, that `TRTSDT` and `TRTEDT` are derived, the function
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_trtdurd.md)
can be used to calculate the Treatment duration (`TRTDURD`).

``` r
adsl <- adsl %>%
  derive_var_trtdurd()
```

### Derive Disposition Variables

#### Disposition Dates (e.g. `EOSDT`)

The functions
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
and
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
can be used to derive a disposition date. First the character
disposition date (`DS.DSSTDTC`) is converted to a numeric date
(`DSSTDT`) calling
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md).
The `DS` dataset is extended by the `DSSTDT` variable because the date
is required by other derivations, e.g., `RANDDT` as well. Then the
relevant disposition date is selected by adjusting the `filter_add`
argument.

To add the End of Study date (`EOSDT`) to the input dataset, a call
could be:

``` r
# Convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )
```

The `ds_ext` dataset:

The `adsl` dataset:

The
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
function allows to impute partial dates as well. If imputation is needed
and missing days are to be imputed to the first of the month and missing
months to the first month of the year, set `highest_imputation = "M"`.

#### Disposition Status (e.g. `EOSSTT`)

The function
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
can be used to derive the End of Study status (`EOSSTT`) based on
`DSCAT` and `DSDECOD` from `DS`. The relevant observations are selected
by adjusting the `filter_add` argument. A function mapping `DSDECOD`
values to `EOSSTT` values can be defined and used in the `new_vars`
argument. The mapping for the call below is

- `"COMPLETED"` if `DSDECOD == "COMPLETED"`
- `NA_character_` if `DSDECOD` is `"SCREEN FAILURE"`
- `"DISCONTINUED"` otherwise

Example function `format_eosstt()`:

``` r
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}
```

The customized mapping function `format_eosstt()` can now be passed to
the main function. For subjects without a disposition event the end of
study status is set to `"ONGOING"` by specifying the `missing_values`
argument.

``` r
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )
```

This call would return the input dataset with the column `EOSSTT` added.

If the derivation must be changed, the user can create his/her own
function to map `DSDECOD` to a suitable `EOSSTT` value.

#### Disposition Reason(s) (e.g. `DCSREAS`, `DCSREASP`)

The main reason for discontinuation is usually stored in `DSDECOD` while
`DSTERM` provides additional details regarding subject’s discontinuation
(e.g., description of `"OTHER"`).

The function
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
can be used to derive a disposition reason (along with the details, if
required) at a specific timepoint. The relevant observations are
selected by adjusting the `filter_add` argument.

To derive the End of Study reason(s) (`DCSREAS` and `DCSREASP`), the
function will map `DCSREAS` as `DSDECOD`, and `DCSREASP` as `DSTERM` if
`DSDECOD` is not `"COMPLETED"`, `"SCREEN FAILURE"`, or `NA`, `NA`
otherwise.

``` r
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD, DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      !(DSDECOD %in% c("SCREEN FAILURE", "COMPLETED", NA))
  )
```

This call would return the input dataset with the column `DCSREAS` and
`DCSREASP` added.

If the derivation must be changed, the user can define that derivation
in the `filter_add` argument of the function to map `DSDECOD` and
`DSTERM` to a suitable `DCSREAS`/`DCSREASP` value.

The call below maps `DCSREAS` and `DCREASP` as follows:

- `DCSREAS` as `DSDECOD` if `DSDECOD` is not `"COMPLETED"` or `NA`, `NA`
  otherwise
- `DCSREASP` as `DSTERM` if `DSDECOD` is equal to `OTHER`, `NA`
  otherwise

``` r
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      DSDECOD %notin% c("SCREEN FAILURE", "COMPLETED", NA)
  ) %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD %in% "OTHER"
  )
```

#### Randomization Date (`RANDDT`)

The function
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
can be used to derive randomization date variable. To map Randomization
Date (`RANDDT`), the call would be:

``` r
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT)
  )
```

This call would return the input dataset with the column `RANDDT` is
added.

### Derive Birth Date and Analysis Age (`BRTHDT`, `AAGE`, `AAGEU`)

The function
[`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_aage.md)
can be used to derive analysis age (`AAGE`) and analysis age unit
(`AAGEU`). The function derives the age based on the birth date
(`BRTHDT`) and a reference date, which is typically the randomization
date (`RANDDT`).

Note that `BRTHDT` must be derived first from `BRTHDTC` using
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md).

``` r
# Derive birth date from BRTHDTC
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "BRTH",
    dtc = BRTHDTC
  )
```

Typically, dates of birth are collected only as years. However, this
data has complete birth dates. The function
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
has the flexibility to do imputation on partial dates. Please see the
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
documentation and examples for more details or check out the [Date and
Time
Imputation](https:/pharmaverse.github.io/admiral/v1.4.1/articles/imputation.md)
vignette.

Now that we have `BRTHDT`, we can use `derive_var_aage()` to derive
`AAGE` and `AAGEU`.

``` r
adsl <- adsl %>%
  derive_vars_aage(
    start_date = BRTHDT,
    end_date = RANDDT
  )
```

This call returns the input dataset with `AAGE` and `AAGEU` added. By
default, the age is calculated in years.

### Derive Death Variables

#### Death Date (`DTHDT`)

The function
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
can be used to derive `DTHDT`. This function allows the user to impute
the date as well. If you have partial dates and are in need of
imputation then please see the
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
documentation and examples for more details or check out the [Date and
Time
Imputation](https:/pharmaverse.github.io/admiral/v1.4.1/articles/imputation.md)
vignette.

``` r
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC
  )
```

#### Cause of Death (`DTHCAUS`)

The cause of death `DTHCAUS` can be derived using the function
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md).

Since the cause of death could be collected/mapped in different domains
(e.g. `DS`, `AE`, `DD`), it is important the user specifies the right
source(s) to derive the cause of death from.

For example, if the date of death is collected in the AE form when the
AE is Fatal, the cause of death would be set to the preferred term
(`AEDECOD`) of that Fatal AE, while if the date of death is collected in
the `DS` form, the cause of death would be set to the disposition term
(`DSTERM`). To achieve this, the
[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md)
objects within
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md)
must be specified and defined such that they fit the study requirement.

An example call to
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md)
would be:

``` r
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS)
  )
```

The function also offers the option to add some traceability variables
(e.g. `DTHDOM` would store the domain where the date of death is
collected, and `DTHSEQ` would store the `xxSEQ` value of that domain).
The traceability variables should be added to the
[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md)
calls and included in the `new_vars` parameter of
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md).

``` r
adsl <- adsl %>%
  select(-DTHCAUS) %>% # Remove it before deriving it again
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = "AE", DTHSEQ = AESEQ),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = "DS", DTHSEQ = DSSEQ),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS, DTHDOM, DTHSEQ)
  )
```

Following the derivation of `DTHCAUS` and related traceability
variables, it is then possible to derive grouping variables such as
death categories (`DTHCGRx`) using standard tidyverse code.

``` r
adsl <- adsl %>%
  mutate(DTHCGR1 = case_when(
    is.na(DTHDOM) ~ NA_character_,
    DTHDOM == "AE" ~ "ADVERSE EVENT",
    str_detect(DTHCAUS, "(PROGRESSIVE DISEASE|DISEASE RELAPSE)") ~ "PROGRESSIVE DISEASE",
    TRUE ~ "OTHER"
  ))
```

#### Duration Relative to Death

The function
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md)
can be used to derive duration relative to death like the Relative Day
of Death (`DTHADY`) or the numbers of days from last dose to death
(`LDDTHELD`).

Example calls:

- Relative Day of Death

``` r
adsl <- adsl %>%
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  )
```

- Elapsed Days from Last Dose to Death

``` r
adsl <- adsl %>%
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )
```

### Derive the Last Date Known Alive (`LSTALVDT`)

Similarly as for the cause of death (`DTHCAUS`), the last known alive
date (`LSTALVDT`) can be derived from multiples sources using
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md).

An example could be (`--DTC` dates are converted to numeric dates
imputing missing day and month to the first):

``` r
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          seq = LBSEQ
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, seq = 0),
      )
    ),
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, seq, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT)
  )
```

Traceability variables can be added by specifying the variables in the
`set_values_to` parameter of the
[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md)
function.

``` r
adsl <- adsl %>%
  select(-LSTALVDT) %>% # Created in the previous call
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AESTDTC"
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AEENDTC"
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          LALVSEQ = LBSEQ,
          LALVDOM = "LB",
          LALVVAR = "LBDTC"
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, LALVSEQ = NA_integer_, LALVDOM = "ADSL", LALVVAR = "TRTEDTM"),
      )
    ),
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, LALVSEQ, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT, LALVSEQ, LALVDOM, LALVVAR)
  )
```

### Derive Groupings and Populations

#### Grouping (e.g. `AGEGR1` or `REGION1`)

Numeric and categorical variables (`AGE`, `RACE`, `COUNTRY`, etc.) may
need to be grouped to perform the required analysis.
[admiral](https://pharmaverse.github.io/admiral/) provides the
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md)
function to create such groups. This function is especially useful if
more than one variable needs to be created for each condition, e.g.,
`AGEGR1` and `AGEGR1N`.

Additionally, one needs to be careful when considering the order of the
conditions in the lookup table. The category is assigned based on the
first match. That means *catch-all* conditions must come after specific
conditions, e.g. `!is.na(COUNTRY)` must come after
`COUNTRY %in% c("CAN", "USA")`.

``` r
# Create lookup tables
agegr1_lookup <- exprs(
  ~condition,           ~AGEGR1,
  AGE < 18,               "<18",
  between(AGE, 18, 64), "18-64",
  AGE > 64,               ">64",
  is.na(AGE),         "Missing"
)

region1_lookup <- exprs(
  ~condition,                          ~REGION1,
  COUNTRY %in% c("CAN", "USA"), "North America",
  !is.na(COUNTRY),          "Rest of the World",
  is.na(COUNTRY),                     "Missing"
)
```

``` r
adsl <- adsl %>%
  derive_vars_cat(
    definition = agegr1_lookup
  ) %>%
  derive_vars_cat(
    definition = region1_lookup
  )
```

Alternatively, you can also solve this task with custom functions:

``` r
format_agegr1 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 64) ~ "18-64",
    var_input > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}
format_region1 <- function(var_input) {
  case_when(
    var_input %in% c("CAN", "USA") ~ "North America",
    !is.na(var_input) ~ "Rest of the World",
    TRUE ~ "Missing"
  )
}

adsl %>%
  mutate(
    AGEGR1 = format_agegr1(AAGE),
    REGION1 = format_region1(COUNTRY)
  )
```

#### Population Flags (e.g. `SAFFL`)

Since the populations flags are mainly company/study specific no
dedicated functions are provided, but in most cases they can easily be
derived using `derive_var_merged_exist_flag`.

An example of an implementation could be:

``` r
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    false_value = "N",
    missing_value = "N",
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  )
```

### Derive Other Variables

The users can add specific code to cover their need for the analysis.

The following functions are helpful for many ADSL derivations:

- [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md) -
  Merge Variables from a Dataset to the Input Dataset
- [`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_exist_flag.md) -
  Merge an Existence Flag
- [`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_summary.md) -
  Merge Summary Variables

See also [Generic
Functions](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

### Add Labels and Attributes

Note that attributes may not be preserved in some cases after processing
with [admiral](https://pharmaverse.github.io/admiral/). The recommended
approach is to apply variable labels and other metadata as a final step
in your data derivation process using packages like:

- [metacore](https://atorus-research.github.io/metacore/): establish a
  common foundation for the use of metadata within an R session.

- [metatools](https://pharmaverse.github.io/metatools/): enable the use
  of metacore objects. Metatools can be used to build datasets or
  enhance columns in existing datasets as well as checking datasets
  against the metadata.

- [xportr](https://atorus-research.github.io/xportr/): functionality to
  associate all metadata information to a local R data frame, perform
  data set level validation checks and convert into a [transport v5
  file(xpt)](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/movefile/n1xbwdre0giahfn11c99yjkpi2yb.htm).

NOTE: Together with [admiral](https://pharmaverse.github.io/admiral/)
these packages comprise an End to End pipeline under the umbrella of the
[pharmaverse](https://github.com/pharmaverse). An example of applying
metadata and perform associated checks can be found at the [pharmaverse
E2E example](https://pharmaverse.github.io/examples/adam/adsl).

## Example Script

| ADaM | Sourcing Command          |
|------|---------------------------|
| ADSL | `use_ad_template("ADSL")` |
