# Creating a BDS Finding ADaM

## Introduction

This article describes creating a BDS finding ADaM. Examples are
currently presented and tested in the context of ADVS. However, the
examples could be applied to other BDS Finding ADaMs such as ADEG, ADLB,
etc. where a single result is captured in an SDTM Finding domain on a
single date and/or time.

**Note**: *All examples assume CDISC SDTM and/or ADaM format as input
unless otherwise specified.*

## Programming Workflow

- [Read in Data](#readdata)
- [Derive/Impute Numeric Date/Time and Analysis Day (`ADT`, `ADTM`,
  `ADY`, `ADTF`, `ATMF`)](#datetime)
- [Assign `PARAMCD`, `PARAM`, `PARAMN`, `PARCAT1`](#paramcd)
- [Derive Results (`AVAL`, `AVALC`)](#aval)
- [Derive Additional Parameters (e.g. `BSA`, `BMI`, or `MAP` for
  `ADVS`)](#derive_param)
- [Derive Timing Variables (e.g. `APHASE`, `AVISIT`,
  `APERIOD`)](#timing)
- [Timing Flag Variables (e.g. `ONTRTFL`)](#timingflag)
- [Assign Reference Range Indicator (`ANRIND`)](#referencerange)
- [Derive Baseline (`BASETYPE`, `ABLFL`, `BASE`, `BASEC`,
  `BNRIND`)](#baseline)
- [Derive Change from Baseline (`CHG`, `PCHG`)](#bchange)
- [Derive Shift (e.g.`SHIFT1`)](#shift)
- [Derive Analysis Ratio (e.g. `R2BASE`)](#analysisratio)
- [Derive Analysis Flags (e.g. `ANL01FL`)](#analysisrec)
- [Assign Treatment (`TRTA`, `TRTP`)](#treatment)
- [Derive Categorization Variables (`AVALCATy`)](#cat)
- [Derive Criterion Variables (`CRITy`, `CRITyFL`,
  `CRITyFN`)](#crit_vars)
- [Derive New Rows](#additional)
- [Assign `ASEQ`](#aseq)
- [Add ADSL variables](#adsl_vars)
- [Add Labels and Attributes](#attributes)

### Read in Data

To start, all data frames needed for the creation of `ADVS` should be
read into the environment. This will be a company specific process. Some
of the data frames needed may be `VS` and `ADSL`.

For example purpose, the CDISC Pilot SDTM and ADaM datasets—which are
included in
[pharmaversesdtm](https://pharmaverse.github.io/pharmaversesdtm/)—are
used.

``` r
library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)
library(tibble)

vs <- pharmaversesdtm::vs
adsl <- admiral::admiral_adsl

vs <- convert_blanks_to_na(vs)
```

At this step, it may be useful to join `ADSL` to your `VS` domain. Only
the `ADSL` variables used for derivations are selected at this step. The
rest of the relevant `ADSL` variables would be added later.

``` r
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

advs <- derive_vars_merged(
  vs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = exprs(STUDYID, USUBJID)
)
```

### Derive/Impute Numeric Date/Time and Analysis Day (`ADT`, `ADTM`, `ADY`, `ADTF`, `ATMF`)

The function
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
can be used to derive `ADT`. This function allows the user to impute the
date as well.

Example calls:

``` r
advs <- derive_vars_dt(advs, new_vars_prefix = "A", dtc = VSDTC)
```

If imputation is needed and the date is to be imputed to the first of
the month, the call would be:

``` r
advs <- derive_vars_dt(
  advs,
  new_vars_prefix = "A",
  dtc = VSDTC,
  highest_imputation = "M"
)
```

Similarly, `ADTM` may be created using the function
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm.md).
Imputation may be done on both the date and time components of `ADTM`.

``` r
# CDISC Pilot data does not contain times and the output of the derivation
# ADTM is not presented.
advs <- derive_vars_dtm(
  advs,
  new_vars_prefix = "A",
  dtc = VSDTC,
  highest_imputation = "M"
)
```

By default, the variable `ADTF` for
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
or `ADTF` and `ATMF` for
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm.md)
will be created and populated with the controlled terminology outlined
in the ADaM IG for date imputations.

See also [Date and Time
Imputation](https:/pharmaverse.github.io/admiral/test_cicd/articles/imputation.md).

Once `ADT` is derived, the function
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dy.md)
can be used to derive `ADY`. This example assumes both `ADT` and
`TRTSDT` exist on the data frame.

``` r
advs <-
  derive_vars_dy(advs, reference_date = TRTSDT, source_vars = exprs(ADT))
```

### Assign `PARAMCD`, `PARAM`, `PARAMN`, `PARCAT1`

To assign parameter level values such as `PARAMCD`, `PARAM`, `PARAMN`,
`PARCAT1`, etc., a lookup can be created to join to the source data.

For example, when creating `ADVS`, a lookup based on the SDTM `--TESTCD`
value may be created:

| `VSTESTCD` | `PARAMCD` | `PARAM`                         | `PARAMN` | `PARCAT1`              | `PARCAT1N` |
|------------|-----------|---------------------------------|----------|------------------------|------------|
| HEIGHT     | HEIGHT    | Height (cm)                     | 1        | Subject Characteristic | 1          |
| WEIGHT     | WEIGHT    | Weight (kg)                     | 2        | Subject Characteristic | 1          |
| DIABP      | DIABP     | Diastolic Blood Pressure (mmHg) | 3        | Vital Sign             | 2          |
| MAP        | MAP       | Mean Arterial Pressure          | 4        | Vital Sign             | 2          |
| PULSE      | PULSE     | Pulse Rate (beats/min)          | 5        | Vital Sign             | 2          |
| SYSBP      | SYSBP     | Systolic Blood Pressure (mmHg)  | 6        | Vital Sign             | 2          |
| TEMP       | TEMP      | Temperature (C)                 | 7        | Vital Sign             | 2          |

This lookup may now be joined to the source data:

At this stage, only `PARAMCD` is required to perform the derivations.
Additional derived parameters may be added, so only `PARAMCD` is joined
to the datasets at this point. All other variables related to `PARAMCD`
(e.g. `PARAM`, `PARCAT1`, …) will be added when all `PARAMCD` are
derived.

``` r
advs <- derive_vars_merged_lookup(
  advs,
  dataset_add = param_lookup,
  new_vars = exprs(PARAMCD),
  by_vars = exprs(VSTESTCD)
)
#> All `VSTESTCD` are mapped.
```

Please note, it may be necessary to include other variables in the join.
For example, perhaps the `PARAMCD` is based on `VSTESTCD` and `VSPOS`,
it may be necessary to expand this lookup or create a separate look up
for `PARAMCD`.

If more than one lookup table, e.g., company parameter mappings and
project parameter mappings, are available,
[`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/consolidate_metadata.md)
can be used to consolidate these into a single lookup table.

Additionally note that each parameter is mapped to only one `PARCAT1`
variable. This is described in section 3.3.4.1 of the ADaM
Implementation Guide, version 1.3: “Any given `PARAM` may be associated
with at-most one level of `PARCATy` (e.g., one level of `PARCAT1` and
one level of `PARCAT2`)”.

### Derive Results (`AVAL`, `AVALC`)

The mapping of `AVAL` and `AVALC` is left to the ADaM programmer. An
example mapping may be:

``` r
advs <- mutate(
  advs,
  AVAL = VSSTRESN
)
```

In this example, as is often the case for ADVS, all `AVAL` values are
numeric without any corresponding non-redundant text value for `AVALC`.
Per recommendation in ADaMIG v1.3 we do not map `AVALC`.

### Derive Additional Parameters (e.g. `BSA`, `BMI` or `MAP` for `ADVS`)

Optionally derive new parameters creating `PARAMCD` and `AVAL`. Note
that only variables specified in the `by_vars` argument will be
populated in the newly created records. This is relevant to the
functions `derive_param_map`, `derive_param_bsa`, `derive_param_bmi`,
and `derive_param_qtc`.

Below is an example of creating `Mean Arterial Pressure` for `ADVS`, see
also Example 3 in section below [Derive New Rows](#additional) for
alternative way of creating new parameters.

``` r
advs <- derive_param_map(
  advs,
  by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
  set_values_to = exprs(PARAMCD = "MAP"),
  get_unit_expr = VSSTRESU,
  filter = VSSTAT != "NOT DONE" | is.na(VSSTAT)
)
```

Likewise, function call below, to create parameter `Body Surface Area`
(BSA) and `Body Mass Index` (BMI) for `ADVS` domain. Note that if height
is collected only once use `constant_by_vars` to specify the
subject-level variable to merge on. Otherwise BSA and BMI are only
calculated for visits where both are collected.

``` r
advs <- derive_param_bsa(
  advs,
  by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
  method = "Mosteller",
  set_values_to = exprs(PARAMCD = "BSA"),
  get_unit_expr = VSSTRESU,
  filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
  constant_by_vars = exprs(USUBJID)
)

advs <- derive_param_bmi(
  advs,
  by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
  set_values_to = exprs(PARAMCD = "BMI"),
  get_unit_expr = VSSTRESU,
  filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
  constant_by_vars = exprs(USUBJID)
)
```

Similarly, for `ADEG`, the parameters `QTCBF` `QTCBS` and `QTCL` can be
created with a function call. See example below for `PARAMCD` = `QTCF`.

``` r
adeg <- tibble::tribble(
  ~USUBJID, ~EGSTRESU, ~PARAMCD, ~AVAL,          ~VISIT,
  "P01",       "msec",     "QT",   350, "CYCLE 1 DAY 1",
  "P01",       "msec",     "QT",   370, "CYCLE 2 DAY 1",
  "P01",       "msec",     "RR",   842, "CYCLE 1 DAY 1",
  "P01",       "msec",     "RR",   710, "CYCLE 2 DAY 1"
)

adeg <- derive_param_qtc(
  adeg,
  by_vars = exprs(USUBJID, VISIT),
  method = "Fridericia",
  set_values_to = exprs(PARAMCD = "QTCFR"),
  get_unit_expr = EGSTRESU
)
```

Similarly, for `ADLB`, the function
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_wbc_abs.md)
can be used to create new parameter for lab differentials converted to
absolute values. See example below:

``` r
adlb <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVAL,                        ~PARAM,          ~VISIT,
  "P01",       "WBC",    33,    "Leukocyte Count (10^9/L)", "CYCLE 1 DAY 1",
  "P01",       "WBC",    38,    "Leukocyte Count (10^9/L)", "CYCLE 2 DAY 1",
  "P01",     "LYMLE",  0.90, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
  "P01",     "LYMLE",  0.70, "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1"
)

derive_param_wbc_abs(
  dataset = adlb,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(
    PARAMCD = "LYMPH",
    PARAM = "Lymphocytes Abs (10^9/L)",
    DTYPE = "CALCULATION"
  ),
  get_unit_expr = extract_unit(PARAM),
  wbc_code = "WBC",
  diff_code = "LYMLE",
  diff_type = "fraction"
)
```

When all `PARAMCD` have been derived and added to the dataset, the other
information from the look-up table (`PARAM`, `PARAMCAT1`,…) should be
added.

``` r
# Derive PARAM and PARAMN
advs <- derive_vars_merged(
  advs,
  dataset_add = select(param_lookup, -VSTESTCD),
  by_vars = exprs(PARAMCD)
)
```

### Derive Timing Variables (e.g. `APHASE`, `AVISIT`, `APERIOD`)

Categorical timing variables are protocol and analysis dependent. Below
is a simple example.

``` r
advs <- advs %>%
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ NA_character_,
      str_detect(VISIT, "UNSCHED") ~ NA_character_,
      str_detect(VISIT, "RETRIEVAL") ~ NA_character_,
      str_detect(VISIT, "AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT)
    ),
    AVISITN = as.numeric(case_when(
      VISIT == "BASELINE" ~ "0",
      str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", ""))
    )),
    ATPT = VSTPT,
    ATPTN = VSTPTNUM
  )

count(advs, VISITNUM, VISIT, AVISITN, AVISIT)
#> # A tibble: 15 × 5
#>    VISITNUM VISIT               AVISITN AVISIT       n
#>       <dbl> <chr>                 <dbl> <chr>    <int>
#>  1      1   SCREENING 1              NA NA         102
#>  2      2   SCREENING 2              NA NA          78
#>  3      3   BASELINE                  0 Baseline    96
#>  4      3.5 AMBUL ECG PLACEMENT      NA NA          65
#>  5      4   WEEK 2                    2 Week 2      96
#>  6      5   WEEK 4                    4 Week 4      80
#>  7      6   AMBUL ECG REMOVAL        NA NA          52
#>  8      7   WEEK 6                    6 Week 6      48
#>  9      8   WEEK 8                    8 Week 8      48
#> 10      9   WEEK 12                  12 Week 12     48
#> 11     10   WEEK 16                  16 Week 16     48
#> 12     11   WEEK 20                  20 Week 20     32
#> 13     12   WEEK 24                  24 Week 24     32
#> 14     13   WEEK 26                  26 Week 26     32
#> 15    201   RETRIEVAL                NA NA          26

count(advs, VSTPTNUM, VSTPT, ATPTN, ATPT)
#> # A tibble: 4 × 5
#>   VSTPTNUM VSTPT                          ATPTN ATPT                           n
#>      <dbl> <chr>                          <dbl> <chr>                      <int>
#> 1      815 AFTER LYING DOWN FOR 5 MINUTES   815 AFTER LYING DOWN FOR 5 MI…   232
#> 2      816 AFTER STANDING FOR 1 MINUTE      816 AFTER STANDING FOR 1 MINU…   232
#> 3      817 AFTER STANDING FOR 3 MINUTES     817 AFTER STANDING FOR 3 MINU…   232
#> 4       NA NA                                NA NA                           187
```

For assigning visits based on time windows and deriving periods,
subperiods, and phase variables see the [“Visit and Period Variables”
vignette](https:/pharmaverse.github.io/admiral/test_cicd/articles/visits_periods.md).

### Timing Flag Variables (e.g. `ONTRTFL`)

In some analyses, it may be necessary to flag an observation as
on-treatment. The admiral function
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_ontrtfl.md)
can be used.

For example, if on-treatment is defined as any observation between
treatment start and treatment end, the flag may be derived as:

``` r
advs <- derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT
)
```

This function returns the original data frame with the column `ONTRTFL`
added. Additionally, this function does have functionality to handle a
window on the `ref_end_date`. For example, if on-treatment is defined as
between treatment start and treatment end plus 60 days, the call would
be:

``` r
advs <- derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  ref_end_window = 60
)
```

In addition, the function does allow you to filter out pre-treatment
observations that occurred on the start date. For example, if
observations with `VSTPT == PRE` should not be considered on-treatment
when the observation date falls between the treatment start and end
date, the user may specify this using the `filter_pre_timepoint`
parameter:

``` r
advs <- derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  filter_pre_timepoint = ATPT == "AFTER LYING DOWN FOR 5 MINUTES"
)
```

Lastly, the function does allow you to create any on-treatment flag
based on the analysis needs. For example, if variable `ONTR01FL` is
needed, showing the on-treatment flag during Period 01, you need to set
`new var = ONTR01FL`. In addition, for Period 01 Start Date and Period
01 End Date, you need `ref_start_date = AP01SDT` and
`ref_end_date = AP01EDT`.

``` r
advs <- derive_var_ontrtfl(
  advs,
  new_var = ONTR01FL,
  start_date = ASTDT,
  end_date = AENDT,
  ref_start_date = AP01SDT,
  ref_end_date = AP01EDT,
  span_period = TRUE
)
```

### Assign Reference Range Indicator (`ANRIND`)

The admiral function
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_anrind.md)
may be used to derive the reference range indicator `ANRIND`.

This function requires the reference range boundaries to exist on the
data frame (`ANRLO`, `ANRHI`) and also accommodates the additional
boundaries `A1LO` and `A1HI`.

The function is called as:

``` r
advs <- derive_var_anrind(advs)
```

### Derive Baseline (`BASETYPE`, `ABLFL`, `BASE`, `BNRIND`)

The `BASETYPE` should be derived using the function
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md).
The parameter `basetypes` of this function requires a named list of
expression detailing how the `BASETYPE` should be assigned. Note, if a
record falls into multiple expressions within the basetypes expression,
a row will be produced for each `BASETYPE`.

``` r
advs <- derive_basetype_records(
  dataset = advs,
  basetypes = exprs(
    "LAST: AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
    "LAST: AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
    "LAST: AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
    "LAST" = is.na(ATPTN)
  )
)

count(advs, ATPT, ATPTN, BASETYPE)
#> # A tibble: 4 × 4
#>   ATPT                           ATPTN BASETYPE                                n
#>   <chr>                          <dbl> <chr>                               <int>
#> 1 AFTER LYING DOWN FOR 5 MINUTES   815 LAST: AFTER LYING DOWN FOR 5 MINUT…   232
#> 2 AFTER STANDING FOR 1 MINUTE      816 LAST: AFTER STANDING FOR 1 MINUTE     232
#> 3 AFTER STANDING FOR 3 MINUTES     817 LAST: AFTER STANDING FOR 3 MINUTES    232
#> 4 NA                                NA LAST                                  187
```

It is important to derive `BASETYPE` first so that it can be utilized in
subsequent derivations. This will be important if the data frame
contains multiple values for `BASETYPE`.

Next, the analysis baseline flag `ABLFL` can be derived using the
[admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md).
For example, if baseline is defined as the last non-missing `AVAL` prior
or on `TRTSDT`, the function call for `ABLFL` would be:

``` r
advs <- restrict_derivation(
  advs,
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = exprs(ADT, ATPTN, VISITNUM),
    new_var = ABLFL,
    mode = "last"
  ),
  filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
)
```

Note: Additional examples of the
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md)
function can be found [above.](#analysisrec)

Lastly, the `BASE`, and `BNRIND` columns can be derived using the
[admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_base.md).
Example calls are:

``` r
advs <- derive_var_base(
  advs,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
  source_var = AVAL,
  new_var = BASE
)

advs <- derive_var_base(
  advs,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
  source_var = ANRIND,
  new_var = BNRIND
)
```

### Derive Change from Baseline (`CHG`, `PCHG`)

Change and percent change from baseline can be derived using the
[admiral](https://pharmaverse.github.io/admiral/) functions
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md)
and
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md).
These functions expect `AVAL` and `BASE` to exist in the data frame. The
`CHG` is simply `AVAL - BASE` and the `PCHG` is
`(AVAL - BASE) / absolute value (BASE) * 100`. Examples calls are:

``` r
advs <- derive_var_chg(advs)

advs <- derive_var_pchg(advs)
```

If the variables should not be derived for all records, e.g., for
post-baseline records only,
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_derivation.md)
can be used.

### Derive Shift (e.g. `SHIFT1`)

Shift variables can be derived using the
[admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md).
This function derives a character shift variable concatenating shift in
values based on a user-defined pairing, e.g., shift from baseline
reference range `BNRIND` to analysis reference range `ANRIND`. Examples
calls are:

``` r
advs <- derive_var_shift(advs,
  new_var = SHIFT1,
  from_var = BNRIND,
  to_var = ANRIND
)
```

If the variables should not be derived for all records, e.g., for
post-baseline records only,
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_derivation.md)
can be used.

### Derive Analysis Ratio (`R2BASE`)

Analysis ratio variables can be derived using the
[admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md).
This function derives a ratio variable based on user-specified pair. For
example, Ratio to Baseline is calculated by `AVAL / BASE` and the
function appends a new variable `R2BASE` to the dataset. This function
can also derive `R2AyHI` and `R2AyLO` values. Examples calls are:

``` r
advs <- derive_var_analysis_ratio(advs,
  numer_var = AVAL,
  denom_var = BASE
)

advs <- derive_var_analysis_ratio(advs,
  numer_var = AVAL,
  denom_var = ANRLO,
  new_var = R01ANRLO
)
```

If the variables should not be derived for all records, e.g., for
post-baseline records only,
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_derivation.md)
can be used.

### Derive Analysis Flags (e.g. `ANL01FL`)

In most finding ADaMs, an analysis flag is derived to identify the
appropriate observation(s) to use for a particular analysis when a
subject has multiple observations within a particular timing period.

In this situation, an analysis flag (e.g. `ANLxxFL`) may be used to
choose the appropriate record for analysis.

This flag may be derived using the
[admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md).
For this example, we will assume we would like to choose the latest and
highest value by `USUBJID`, `PARAMCD`, `AVISIT`, and `ATPT`.

``` r
advs <- restrict_derivation(
  advs,
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD, AVISIT),
    order = exprs(ADT, ATPTN, AVAL),
    new_var = ANL01FL,
    mode = "last"
  ),
  filter = !is.na(AVISITN)
)
```

Another common example would be flagging the worst value for a subject,
parameter, and visit. For this example, we will assume we have 3
`PARAMCD` values (`SYSBP`, `DIABP`, and `RESP`). We will also assume
high is worst for `SYSBP` and `DIABP` and low is worst for `RESP`.

``` r
advs <- slice_derivation(
  advs,
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD, AVISIT),
    order = exprs(ADT, ATPTN),
    new_var = WORSTFL,
    mode = "first"
  ),
  derivation_slice(
    filter = PARAMCD %in% c("SYSBP", "DIABP") & (!is.na(AVISIT) & !is.na(AVAL))
  ),
  derivation_slice(
    filter = PARAMCD %in% "PULSE" & (!is.na(AVISIT) & !is.na(AVAL)),
    args = params(mode = "last")
  )
) %>%
  arrange(STUDYID, USUBJID, BASETYPE, PARAMCD, AVISIT)
```

### Assign Treatment (`TRTA`, `TRTP`)

`TRTA` and `TRTP` must match at least one value of the character
treatment variables in ADSL (e.g., `TRTxxA`/`TRTxxP`,
`TRTSEQA`/`TRTSEQP`, `TRxxAGy`/`TRxxPGy`).

An example of a simple implementation for a study without periods could
be:

``` r
advs <- mutate(advs, TRTP = TRT01P, TRTA = TRT01A)

count(advs, TRTP, TRTA, TRT01P, TRT01A)
#> # A tibble: 2 × 5
#>   TRTP                TRTA                TRT01P              TRT01A           n
#>   <chr>               <chr>               <chr>               <chr>        <int>
#> 1 Placebo             Placebo             Placebo             Placebo        640
#> 2 Xanomeline Low Dose Xanomeline Low Dose Xanomeline Low Dose Xanomeline …   243
```

For studies with periods see the [“Visit and Period Variables”
vignette](https:/pharmaverse.github.io/admiral/test_cicd/articles/visits_periods.html#treatment_bds).

### Derive Categorization Variables (`AVALCATy`)

We can use the
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_cat.md)
function to derive the categorization variables.

``` r
avalcat_lookup <- exprs(
  ~PARAMCD,  ~condition,   ~AVALCAT1, ~AVALCA1N,
  "HEIGHT",  AVAL > 140,   ">140 cm",         1,
  "HEIGHT", AVAL <= 140, "<= 140 cm",         2
)
advs <- advs %>%
  derive_vars_cat(
    definition = avalcat_lookup,
    by_vars = exprs(PARAMCD)
  )
```

### Derive Criterion Variables (`CRITy`, `CRITyFL`, `CRITyFN`)

For deriving criterion variables (`CRITy`, `CRITyFL`, `CRITyFN`)
[admiral](https://pharmaverse.github.io/admiral/) provides
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_crit_flag.md).
It ensures that they are derived in an ADaM-compliant way (see
documentation of the function for details).

In most cases the criterion depends on the parameter. The higher order
functions
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_derivation.md)
and
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/slice_derivation.md)
are useful in this case. In the following example the criterion flags
for systolic and diastolic blood pressure from the ADaM IG are derived.

The first criterion is based on `AVAL` and is derived for systolic and
diastolic blood pressure.
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/slice_derivation.md)
us used to specify the condition and description of the criterion
depending on the parameter.

``` r
advs <- advs %>%
  slice_derivation(
    derivation = derive_vars_crit_flag,
    args = params(
      values_yn = TRUE,
      create_numeric_flag = TRUE
    ),
    derivation_slice(
      filter = PARAMCD == "SYSBP",
      args = params(
        condition = AVAL > 160,
        description = "Systolic Pressure > 160"
      )
    ),
    derivation_slice(
      filter = PARAMCD == "DIABP",
      args = params(
        condition = AVAL > 95,
        description = "Diastolic Pressure > 95"
      )
    )
  )
```

The second criterion is based on `AVAL` and `CHG` and is derived for
systolic blood pressure only. Thus
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_derivation.md)
is used.

``` r
advs <- advs %>%
  restrict_derivation(
    derivation = derive_vars_crit_flag,
    args = params(
      condition = AVAL > 160 & CHG > 10,
      description = "Systolic Pressure > 160 and Change from Baseline in Systolic Pressure > 10",
      crit_nr = 2,
      values_yn = TRUE,
      create_numeric_flag = TRUE
    ),
    filter = PARAMCD == "SYSBP"
  )
```

### Derive New Rows

When deriving new rows for a data frame, it is essential the programmer
takes time to insert this derivation in the correct location of the
code. The location will vary depending on what previous computations
should be retained on the new record and what computations must be done
with the new records.

#### Example 1 (Creating a New Record):

To add a new record based on the selection of a certain criterion
(e.g. minimum, maximum)
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_records.md)
can be used. The new records include all variables of the selected
records.

##### Adding a New Record for the Last Value

For each subject and Vital Signs parameter, add a record holding last
valid observation before end of treatment. Set `AVISIT` to
`"End of Treatment"` and assign a unique `AVISITN` value.

``` r
advs_ex1 <- advs %>%
  derive_extreme_records(
    dataset_add = advs,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    order = exprs(ADT, AVISITN, ATPTN, AVAL),
    mode = "last",
    filter_add = (4 < AVISITN & AVISITN <= 12 & ANL01FL == "Y"),
    set_values_to = exprs(
      AVISIT = "End of Treatment",
      AVISITN = 99,
      DTYPE = "LOV"
    )
  )
```

##### Adding a New Record for the Minimum Value

For each subject and Vital Signs parameter, add a record holding the
minimum value before end of treatment. If the minimum is attained by
multiple observations the first one is selected. Set `AVISIT` to
`"Minimum on Treatment"` and assign a unique `AVISITN` value.

``` r
advs_ex1 <- advs %>%
  derive_extreme_records(
    dataset_add = advs,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    order = exprs(AVAL, ADT, AVISITN, ATPTN),
    mode = "first",
    filter_add = (4 < AVISITN & AVISITN <= 12 & ANL01FL == "Y" & !is.na(AVAL)),
    set_values_to = exprs(
      AVISIT = "Minimum on Treatment",
      AVISITN = 98,
      DTYPE = "MINIMUM"
    )
  )
```

#### Example 2 (Deriving a Summary Record)

For adding new records based on aggregating records
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)
can be used. For the new records only the variables specified by
`by_vars` and `set_values_to` are populated.

For each subject, Vital Signs parameter, visit, and date add a record
holding the average value for observations on that date. Set `DTYPE` to
`AVERAGE`.

``` r
advs_ex2 <- derive_summary_records(
  advs,
  dataset_add = advs,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD, VISITNUM, ADT),
  set_values_to = exprs(
    AVAL = mean(AVAL, na.rm = TRUE),
    DTYPE = "AVERAGE"
  )
)
```

#### Example 3 (Deriving a New `PARAMCD`)

Use function
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md)
to create a new `PARAMCD`. Note that only variables specified in the
`by_vars` argument will be populated in the newly created records.

Below is an example of creating `Mean Arterial Pressure`
(`PARAMCD = MAP2`) with an alternative formula.

``` r
advs_ex3 <- derive_param_computed(
  advs,
  by_vars = exprs(USUBJID, VISIT, ATPT),
  parameters = c("SYSBP", "DIABP"),
  set_values_to = exprs(
    AVAL = (AVAL.SYSBP - AVAL.DIABP) / 3 + AVAL.DIABP,
    PARAMCD = "MAP2",
    PARAM = "Mean Arterial Pressure 2 (mmHg)"
  )
)
```

### Assign `ASEQ`

The [admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_obs_number.md)
can be used to derive `ASEQ`. An example call is:

``` r
advs <- derive_var_obs_number(
  advs,
  new_var = ASEQ,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN),
  check_type = "error"
)
```

### Add ADSL variables

If needed, the other `ADSL` variables can now be added. List of ADSL
variables already merged held in vector `adsl_vars`

``` r
advs <- advs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )
```

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

## Example Scripts

| ADaM | Sourcing Command          |
|------|---------------------------|
| ADEG | `use_ad_template("ADEG")` |
| ADVS | `use_ad_template("ADVS")` |
| ADLB | `use_ad_template("ADLB")` |
