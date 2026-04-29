# Visit and Period Variables

## Introduction

The derivation of visit variables like `AVISIT`, `AVISITN`, `AWLO`,
`AWHI`, … or period, subperiod, or phase variables like `APERIOD`,
`TRT01A`, `TRT02A`, `ASPER`, `PHSDTM`, `PHEDTM`, … is highly
study-specific. Therefore
[admiral](https://pharmaverse.github.io/admiral/) cannot provide
functions which derive these variables. However, for common scenarios
like visit assignments based on time windows or deriving BDS period
variables from ADSL period variables, functions are provided which
support those derivations.

### Required Packages

The examples of this vignette require the following packages.

``` r
library(admiral)
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
```

## Visit variables (`AVISIT`, `AVISITN`, `AWLO`, `AWHI`, …)

The most common ways of deriving `AVISIT` and `AVISITN` are:

- The variables are set to the collected visits (`VISIT` and
  `VISITNUM`).
- The variables are set based on time windows.

The former can be achieved simply by calling
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html), like in
the vignettes and the template scripts.

For the latter a (study-specific) reference dataset needs to be created
which provides for each visit the start and end day (`AWLO` and `AWHI`)
and the values of other visit related variables (`AVISIT`, `AVISITN`,
`AWTARGET`, …).

``` r
windows <- tribble(
  ~AVISIT,    ~AWLO, ~AWHI, ~AVISITN, ~AWTARGET,
  "BASELINE",   -30,     1,        0,         1,
  "WEEK 1",       2,     7,        1,         5,
  "WEEK 2",       8,    15,        2,        11,
  "WEEK 3",      16,    22,        3,        19,
  "WEEK 4",      23,    30,        4,        26
)
```

Then the visits can be assigned based on the analysis day (`ADY`) by
calling
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md):

``` r
adbds <- tribble(
  ~USUBJID, ~ADY,
  "1",       -33,
  "1",        -2,
  "1",         3,
  "1",        24,
  "2",        NA,
)

adbds1 <- adbds %>%
  derive_vars_joined(
    dataset_add = windows,
    filter_join = AWLO <= ADY & ADY <= AWHI,
    join_type = "all",
  )
```

## Period, Subperiod, and Phase Variables

If periods, subperiods, or phases are used, in the simpler of ADaM
applications the corresponding variables have to be consistent across
all datasets. This can be achieved by defining the periods, subperiods,
or phases once and then using this definition for all datasets. The
definition can be stored in ADSL or in a separate dataset. In this
vignette’s examples, this separate dataset is called period/phase
reference dataset depending what data it contains.

Note that periods, subperiods, or phases *can* be defined differently
across datasets (for instance, they may be different across safety and
efficacy analyses) in which case the start/stop dates should be defined
in the individual datasets, instead of in ADSL. However, this vignette
will not cover this scenario.

### A Note on Study Specific Code

The sections below showcase the available tools in
[admiral](https://pharmaverse.github.io/admiral/) to work with period,
subperiod and phase variables. However, at some point study specific
code will always be required. There are two options:

- Study specific code is used to first derive the variables `PxxSwSDT`
  and `PxxSwEDT` in ADSL. Then
  [`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md)
  and
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md)
  can be used to derive period/subperiod variables like `ASPER` or
  `ASPRSDT` in BDS and OCCDS datasets. See an example dataset
  [here](#adsl_example).

- Study specific code is used to derive a dataset with one observation
  per patient, period, and subperiod (see [period reference
  dataset](#reference)). Then
  [`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md)
  can be used to derive `PxxSwSDT` and `PxxSwEDT` in ADSL and
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md)
  can be used to derive period/subperiod variables like `ASPER` or
  `ASPRSDT` in BDS and OCCDS datasets.

It depends on the specific definition of the periods/subperiods which
option works best. If the definition is based on other ADSL variables,
the first option would work best. If the definition is based on
vertically structured data like exposure data (EX dataset), the second
option should be used. This vignette contains examples for both
workflows.

### The Period/Phase Reference Dataset

The [admiral](https://pharmaverse.github.io/admiral/) functions expect
separate reference datasets for periods, subperiods, and phases. For
periods the numeric variable `APERIOD` is expected, for subperiods the
numeric variables `APERIOD` and `ASPER`, and for phases the numeric
variable `APHASEN`.

The period/phase reference dataset should be created according to the
design and ADaM needs of the study in question. It should contain one
observation per subject and period, subperiod, or phase. See the [next
section](#example_period) for an example dataset.

### Example Creation of the Period/Phase Reference Dataset

Consider a simple crossover study with the following design:

![Dummy Study Design
Image](https://raw.githubusercontent.com/pharmaverse/admiral/main/inst/visit_periods/dummy_study_design.png)

Given this design, an option could be to split the study into two
periods:

- Period 1 (Day 1 to Week 16);
- Period 2 (Week 16 to End of Study).

Alternatively (or additionally) one could split into two phases:

- Phase 1 - Treatment (Day 1 to 28 days after Week 28);
- Phase 2 - Follow-up (29 days after Week 28 to End of Study).

Below, we present two example workflows: one where where a period
reference dataset is created from the exposure dataset EX, and the other
where a phase reference dataset is created using ADSL variables.

#### Creating the Period/Phase Reference Dataset from EX

Below we create a period reference dataset starting from the exposure
dataset EX. Consider the following exposure dataset:

``` r
ex <- tribble(
  ~STUDYID, ~USUBJID,  ~VISIT,    ~EXTRT,   ~EXSTDTC,
  "xyz",    "1",       "Day 1",   "Drug X", "2022-01-02",
  "xyz",    "1",       "Week 4",  "Drug X", "2022-02-05",
  "xyz",    "1",       "Week 8",  "Drug X", "2022-03-01",
  "xyz",    "1",       "Week 12", "Drug X", "2022-04-03",
  "xyz",    "1",       "Week 16", "Drug Y", "2022-05-03",
  "xyz",    "1",       "Week 20", "Drug Y", "2022-06-02",
  "xyz",    "1",       "Week 24", "Drug Y", "2022-07-01",
  "xyz",    "1",       "Week 28", "Drug Y", "2022-08-04",
  "xyz",    "2",       "Day 1",   "Drug Y", "2023-10-20",
  "xyz",    "2",       "Week 4",  "Drug Y", "2023-11-21",
  "xyz",    "2",       "Week 8",  "Drug Y", "2023-12-19",
  "xyz",    "2",       "Week 12", "Drug Y", "2024-01-19",
  "xyz",    "2",       "Week 16", "Drug X", "2024-02-20",
  "xyz",    "2",       "Week 20", "Drug X", "2024-03-17",
  "xyz",    "2",       "Week 24", "Drug X", "2024-04-22",
  "xyz",    "2",       "Week 28", "Drug X", "2024-05-21"
)
```

Then to create a period reference dataset, code like this (or similar)
would suffice:

``` r
period_ref <- ex %>%
  # Select visits marking the start of each period
  filter(VISIT %in% c("Day 1", "Week 16")) %>%
  # Create APERIOD, APERSDT, TRTA based on SDTM counterparts
  mutate(
    APERIOD = case_when(
      VISIT == "Day 1" ~ 1,
      VISIT == "Week 16" ~ 2
    ),
    TRTA = EXTRT,
    APERSDT = convert_dtc_to_dt(EXSTDTC)
  ) %>%
  # Create APEREDT based on start date of next period
  arrange(USUBJID, APERSDT) %>%
  group_by(USUBJID) %>%
  mutate(
    APEREDT = lead(APERSDT) - 1 # one day before start of next period
  ) %>%
  # Tidy up
  ungroup() %>%
  select(-starts_with("EX"), -VISIT)
```

The workflow above populates the Period End Date `APEREDT` for all
periods except the last one. The value of the last period could then be
populated with the End of Study Date (`EOSDT`) from ADSL, to obtain:

``` r
adsl <- tribble(
  ~STUDYID,  ~USUBJID,  ~TRTSDT,           ~TRTEDT,           ~EOSDT,
  "xyz",     "1",       ymd("2022-01-02"), ymd("2022-08-04"), ymd("2022-09-10"),
  "xyz",     "2",       ymd("2023-10-20"), ymd("2024-05-21"), ymd("2024-06-30")
)

period_ref <- period_ref %>%
  left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
  mutate(APEREDT = case_when(
    APERIOD == "1" ~ APEREDT,
    APERIOD == "2" ~ EOSDT
  )) %>%
  select(-EOSDT, -TRTSDT, -TRTEDT)
```

#### Creating the Period/Phase Reference Dataset from ADSL

If the treatment start and end dates are already included in ADSL, we
can derive the phase variables directly in ADSL and create a phase
reference dataset by employing
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md).
Here is an example command to achieve this goal:

``` r
adsl1 <- adsl %>%
  mutate(
    PH1SDT = TRTSDT,
    PH1EDT = TRTEDT + 28,
    APHASE1 = "TREATMENT",
    PH2SDT = TRTEDT + 29,
    PH2EDT = EOSDT,
    APHASE2 = "FUP"
  )

phase_ref <- create_period_dataset(
  adsl1,
  new_vars = exprs(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
)
```

### Creating ADSL Period, Subperiod, or Phase Variables

If a period/phase reference dataset is available, the ADSL variables for
periods, subperiods, or phases can be created from this dataset by
calling
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md).

For example the period reference dataset from the previous section can
be used to add the period variables (`APxxSDT`, `APxxEDT`) to ADSL:

``` r
adsl2 <- derive_vars_period(
  adsl,
  dataset_ref = period_ref,
  new_vars = exprs(APxxSDT = APERSDT, APxxEDT = APEREDT)
)
```

### Creating BDS and OCCDS Period, Subperiod, or Phase Variables

If a period/phase reference dataset is available, BDS and OCCDS
variables for periods, subperiods, or phases can be created by calling
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md).

For example the variables `APHASEN`, `PHSDT`, `PHEDT`, `APHASE` can be
derived from the phase reference dataset defined above.

``` r
adae <- tribble(
  ~STUDYID, ~USUBJID, ~ASTDT,
  "xyz",    "1",      "2022-01-31",
  "xyz",    "1",      "2022-05-02",
  "xyz",    "1",      "2022-09-03",
  "xyz",    "1",      "2022-09-09",
  "xyz",    "2",      "2023-12-25",
  "xyz",    "2",      "2024-06-19",
) %>%
  mutate(ASTDT = ymd(ASTDT))

adae1 <- adae %>%
  derive_vars_joined(
    dataset_add = phase_ref,
    by_vars = exprs(STUDYID, USUBJID),
    filter_join = PHSDT <= ASTDT & ASTDT <= PHEDT,
    join_type = "all"
  )
```

## Treatment Variables (`TRTxxP`, `TRTxxA`, `TRTP`, `TRTA`, …)

In studies with multiple periods the treatment can differ by period,
e.g. for a crossover trial - see the [previous section](#example_period)
for an example design showcasing this. CDISC defines variables for
planned and actual treatments in ADSL (`TRTxxP`, `TRTxxA`, `TRxxPGy`,
`TRxxAGy`, …) and corresponding variables in BDS and OCCDS datasets
(`TRTP`, `TRTA`, `TRTPGy`, `TRTAGy`, …). They can be derived in the same
way (and same step) as the period, subperiod, and phase variables.

### Creating ADSL Treatment Variables

If the treatment information is included in the period/phase reference
dataset, the treatment ADSL variables can be created by calling
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md).
This is showcased below using the period reference dataset from previous
sections.

``` r
adsl <- derive_vars_period(
  adsl,
  dataset_ref = period_ref,
  new_vars = exprs(
    APxxSDT = APERSDT,
    APxxEDT = APEREDT,
    TRTxxA = TRTA
  )
)
```

### Creating BDS and OCCDS Treatment Variables

If a period/phase reference dataset is available, BDS and OCCDS
variables for treatment can be created by calling
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md).

For example the variables `APERIOD` and `TRTA` can be derived from the
period reference dataset defined above.

``` r
adae <- tribble(
  ~STUDYID, ~USUBJID, ~ASTDT,
  "xyz",    "1",      "2022-01-31",
  "xyz",    "1",      "2022-05-02",
  "xyz",    "1",      "2022-08-24",
  "xyz",    "1",      "2022-09-09",
  "xyz",    "2",      "2023-12-25",
  "xyz",    "2",      "2024-06-07",
) %>%
  mutate(ASTDT = ymd(ASTDT))

adae2 <- adae %>%
  derive_vars_joined(
    dataset_add = period_ref,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(APERIOD, TRTA),
    join_vars = exprs(APERSDT, APEREDT),
    join_type = "all",
    filter_join = APERSDT <= ASTDT & ASTDT <= APEREDT
  )
```

If no period/phase reference dataset is available but period/phase
variables are in ADSL, then the former can again be created from ADSL by
calling
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md),
as was showcased [here](#from_adsl).

This time, when calling
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md)
we just need to make sure we include `TRTA = TRTxxA` as part of the
`new_vars` argument to create the treatment variables as well.

``` r
period_ref1 <- adsl %>%
  create_period_dataset(
    new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA)
  )
```
