# Creating a BDS Exposure ADaM

## Introduction

This article describes creating an Exposure ADaM using the BDS
structure. Examples are currently presented using an underlying `EX`
domain where the `EX` domain represents data as collected on the CRF and
the `ADEX` ADaM is output. However, the examples can be applied to
situations where an `EC` domain is used as input and/or `ADEC` or
another exposure ADaM is created.

There are many different approaches to modeling exposure data. This
vignette gives examples of creating `PARAMCD` and `AVAL` combinations
using exposure data. This vignette is not meant to be a guide or
standard for the structure of exposure analysis datasets.

**Note**: *All examples assume CDISC SDTM and/or ADaM format as input
unless otherwise specified.*

## Programming Workflow

- [Read in Data](#readdata)
- [Derive/Impute Numeric Date/Time and Analysis Day (`ADT`, `ADTM`,
  `ADY`, `ADTF`, `ATMF`)](#datetime)
- [Compute Duration for a Record](#duration)
- [Create 1:1 Mapping Records](#onetoone)
- [Create Summary Records](#summaryrec)
- [Assign `PARAMCD`, `PARAMN`, etc. from Reference Tables](#paramcd)
- [Derive Categorization Variables (`AVALCATy`)](#cat)
- [Assign `ASEQ`](#aseq)
- [Add ADSL variables](#adsl_vars)
- [Add Labels and Attributes](#attributes)

### Read in Data

To start, all data frames needed for the creation of `ADEX` should be
read into the environment. This will be a company specific process. Some
of the data frames needed may be `EX` and `ADSL`.

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

ex <- pharmaversesdtm::ex
adsl <- admiral::admiral_adsl

ex <- convert_blanks_to_na(ex)
```

At this step, it may be useful to join `ADSL` to your `EX` domain as
well. Only the `ADSL` variables used for derivations are selected at
this step. The rest of the relevant `ADSL` variables would be added
later.

``` r
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRTEDT, TRTEDTM)

adex <- derive_vars_merged(
  ex,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)
```

The CDISC pilot `EX` domain data does not contain a dose adjustment flag
or the planned dose information. For demonstration purposes, this will
be added to the data.

``` r
adex <- adex %>%
  mutate(
    EXADJ = case_when(
      USUBJID == "01-701-1028" & VISIT %in% c("WEEK 2") ~ "ADVERSE EVENT",
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 2", "WEEK 24") ~ "MEDICATION ERROR",
      TRUE ~ NA_character_
    ),
    EXDOSE = case_when(
      USUBJID == "01-701-1028" & VISIT %in% c("WEEK 2") ~ 0,
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 2", "WEEK 24") ~ 0,
      TRUE ~ EXDOSE
    )
  ) %>%
  mutate(EXPLDOS = if_else(EXTRT == "PLACEBO", 0, 54))

distinct(adex, EXTRT, EXPLDOS)
#> # A tibble: 2 × 2
#>   EXTRT      EXPLDOS
#>   <chr>        <dbl>
#> 1 PLACEBO          0
#> 2 XANOMELINE      54
count(adex, EXADJ)
#> # A tibble: 1 × 2
#>   EXADJ     n
#>   <chr> <int>
#> 1 NA       13
```

### Derive/Impute Numeric Date/Time and Analysis Day (`ADT`, `ADTM`, `ADY`, `ADTF`, `ATMF`)

The function
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
can be used to derive `ADT`. This function allows the user to impute the
date as well.

Example calls:

``` r
adex <- derive_vars_dt(adex, new_vars_prefix = "AST", dtc = EXSTDTC)
adex <- derive_vars_dt(adex, new_vars_prefix = "AEN", dtc = EXENDTC)
```

The next examples demonstrates the datetime imputation features
available in the
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm.md)
function, where the time is imputed as “00:00:00”:

``` r
adex <- derive_vars_dtm(
  adex,
  dtc = EXSTDTC,
  highest_imputation = "M",
  new_vars_prefix = "AST"
)

adex <- derive_vars_dtm(
  adex,
  dtc = EXENDTC,
  highest_imputation = "M",
  date_imputation = "last",
  new_vars_prefix = "AEN"
)
```

The example above imputes the start date to the first first day of the
month and imputes the end date to the last day of the month.

Please see the [Date and Time
Imputation](https:/pharmaverse.github.io/admiral/v1.4.1/articles/imputation.Rmd)
for additional examples on calculating and imputing analysis dates.

Next, the analysis study days can be derived:

``` r
adex <-
  derive_vars_dy(adex,
    reference_date = TRTSDT,
    source_vars = exprs(ASTDT, AENDT)
  )
```

### Compute duration for a record

To compute the duration of treatment or exposure for a record, the
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md)
function can be used.

``` r
adex <- adex %>%
  derive_vars_duration(
    new_var = EXDURD,
    start_date = ASTDT,
    end_date = AENDT
  )
```

The units of the calculated duration can also be changed. In this
example, the duration is output as years:

``` r
adex <- adex %>%
  derive_vars_duration(
    new_var = EXDURDY,
    out_unit = "years",
    start_date = ASTDT,
    end_date = AENDT
  )
```

Please refer to the
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md)
documentation for detailed information on the input parameters.

It may be necessary to calculate additional intermediate values. For
example, the cumulative doses received and cumulative planned doses may
be calculated as:

``` r
adex <- adex %>%
  mutate(
    DOSEO = EXDOSE * EXDURD,
    PDOSEO = EXPLDOS * EXDURD
  )
```

It may be of additional interest to turn a single record containing
dosing summary information into a set of multiple single records, each
representing a single dose over the interval specified by the summary
record. This is another approach to deriving a total dose parameter when
`EXDOSFRQ != ONCE`.

The function
[`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_single_dose_dataset.md)
can be used to expand a record containing a start date, an end date, and
a dosing frequency to a corresponding set of records each representing
one dose (i.e. `EXDOSFRQ == "ONCE"`).

``` r
single_dose <- adex %>%
  filter(USUBJID == "01-701-1015" & EXSTDY == 1) %>%
  create_single_dose_dataset(keep_source_vars = exprs(USUBJID, EXDOSE, EXPLDOS, EXDOSFRQ, ASTDT, AENDT))
```

### Create 1:1 mapping records

The first set of exposure records to create will be records mapped 1:1
to an existing collected exposure record in SDTM. For these records, the
`AVAL` or `AVALC` would be calculated using columns that exist on the
data and no summarizing of records would be necessary.

These records may be used for input into summary records or be used
individually for summarization in outputs. Some examples may be exposure
duration, dose administered, dose adjusted, etc. based on one exposure
record in SDTM.

These records can be derived using simple
[`dplyr::mutate`](https://dplyr.tidyverse.org/reference/mutate.html)
assignments and then combined:

``` r
adex_durd <- adex %>%
  mutate(
    PARAMCD = "DURD",
    AVAL = EXDURD
  )

adex_dose <- adex %>%
  mutate(
    PARAMCD = "DOSE",
    AVAL = DOSEO
  )

adex_pldos <- adex %>%
  mutate(
    PARAMCD = "PLDOSE",
    AVAL = PDOSEO
  )

adex_adj <- adex %>%
  mutate(
    PARAMCD = "ADJ",
    AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)
  )

adex_adjae <- adex %>%
  mutate(
    PARAMCD = "ADJAE",
    AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_)
  )

adex <- bind_rows(
  adex_durd,
  adex_dose,
  adex_pldos,
  adex_adj,
  adex_adjae
) %>%
  mutate(PARCAT1 = "INDIVIDUAL")

count(adex, PARAMCD)
#> # A tibble: 5 × 2
#>   PARAMCD     n
#>   <chr>   <int>
#> 1 ADJ        13
#> 2 ADJAE      13
#> 3 DOSE       13
#> 4 DURD       13
#> 5 PLDOSE     13
```

### Create Summary Records

Exposure is commonly analyzed by a timing interval (e.g. `APHASE`,
`APERIOD`, `AVISIT`, etc.). For these types of calculations, the
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exposure.md)
function may be used. In addition to creating a summarized `AVAL`, the
function will also compute minimum and maximum dates for the record.

For example, to calculate the total dose by subject and treatment,

``` r
adex <- derive_param_exposure(
  adex,
  dataset_add = adex,
  by_vars = c(get_admiral_option("subject_keys"), adsl_vars),
  input_code = "DOSE",
  set_values_to = exprs(
    PARAMCD = "TDOSE",
    PARCAT1 = "OVERALL",
    AVAL = sum(AVAL, na.rm = TRUE)
  )
)
```

A record with `PARAMCD == "TDOSE"` is created with `PARCAT1` set to
`"OVERALL"` using the records in `ADEX` where `PARAMCD == "DOSE"` by
summing `AVAL`. In addition, the `ASTDT`, and `AENDT` are created as the
minimum and maximum date/times associated with each `by_vars` grouping.
Note that, in addition to `PARAMCD`, `PARCAT1`, `AVAL`, `ASTDT` and
`AENDT`, only those variables specified in the `by_vars` argument will
be populated in the new records.

Multiple parameters (records) may be created at one time using the
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md)
function:

``` r
adex <- adex %>%
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = exprs(
          PARAMCD = "TDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TPDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "PLDOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TDURD",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DURD"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJ",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJ"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJAE",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJAE"
      )
    ),
    dataset_add = adex,
    by_vars = c(get_admiral_option("subject_keys"), adsl_vars)
  )

count(adex, PARAMCD, PARCAT1)
#> # A tibble: 10 × 3
#>    PARAMCD PARCAT1        n
#>    <chr>   <chr>      <int>
#>  1 ADJ     INDIVIDUAL    13
#>  2 ADJAE   INDIVIDUAL    13
#>  3 DOSE    INDIVIDUAL    13
#>  4 DURD    INDIVIDUAL    13
#>  5 PLDOSE  INDIVIDUAL    13
#>  6 TADJ    OVERALL        6
#>  7 TADJAE  OVERALL        6
#>  8 TDOSE   OVERALL        6
#>  9 TDURD   OVERALL        6
#> 10 TPDOSE  OVERALL        6
```

Dose intensity can be calculated using the function
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_doseint.md).
The planned dose and administered dose are passed into the function and
a new record is created with the dose intensity calculation. Again, only
those variables specified in the `by_vars` argument will be populated in
this new record.

``` r
adex <- adex %>%
  derive_param_doseint(
    by_vars = c(get_admiral_option("subject_keys"), adsl_vars),
    set_values_to = exprs(PARAMCD = "TNDOSINT"),
    tadm_code = "TDOSE",
    tpadm_code = "TPDOSE"
  )
```

The default calculation for dose intensity is: Administered Doses /
Planned Doses \* 100.

Please see the
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_doseint.md)
documentation to see how planned doses of 0 or `NA` are handled.

### Assign `PARAMCD`, `PARAMN`, etc. from Reference tables

To assign parameter level values such as `PARAM`, `PARAMN`, `PARCAT1`,
etc., a lookup can be created to join to the source data.

For example, when creating `ADEX`, a lookup based on the ADaM `PARAMCD`
value may be created:

| PARAMCD  | PARAM                                                      | PARAMN |
|----------|------------------------------------------------------------|--------|
| DURD     | Study drug duration during constant dosing interval (days) | 1      |
| DOSE     | Dose administered during constant dosing interval (mg)     | 2      |
| PLDOSE   | Planned dose during constant dosing interval (mg)          | 3      |
| ADJ      | Dose adjusted during constant dosing interval              | 4      |
| ADJAE    | Dose adjusted due to AE during constant dosing interval    | 5      |
| TDURD    | Overall duration (days)                                    | 6      |
| TDOSE    | Total dose administered (mg)                               | 7      |
| TPDOSE   | Total planned dose (mg)                                    | 9      |
| TADJ     | Dose adjusted during study                                 | 10     |
| TADJAE   | Dose adjusted during study due to AE                       | 11     |
| TNDOSINT | Overall dose intensity (%)                                 | 12     |

``` r
adex <- derive_vars_merged(
  adex,
  dataset_add = param_lookup,
  by_vars = exprs(PARAMCD)
)

count(adex, PARAMCD, PARAM, PARAMN)
#> # A tibble: 11 × 4
#>    PARAMCD  PARAM                                                   PARAMN     n
#>    <chr>    <chr>                                                    <dbl> <int>
#>  1 ADJ      Dose adjusted during constant dosing interval                4    13
#>  2 ADJAE    Dose adjusted  due to AE during constant dosing interv…      5    13
#>  3 DOSE     Dose administered during constant dosing interval (mg)       2    13
#>  4 DURD     Study drug duration during constant dosing interval (d…      1    13
#>  5 PLDOSE   Planned dose during constant dosing interval (mg)            3    13
#>  6 TADJ     Dose adjusted during study                                  10     6
#>  7 TADJAE   Dose adjusted during study due to AE                        11     6
#>  8 TDOSE    Total dose administered (mg)                                 7     6
#>  9 TDURD    Overall duration (days)                                      6     6
#> 10 TNDOSINT Overall dose intensity (%)                                  12     6
#> 11 TPDOSE   Total planned dose (mg)                                      9     6
```

Please note, this is an example only and additional columns may be
needed for the join depending on your lookup/metadata table.

### Derive Categorization Variables (`AVALCATy`)

We can use the
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md)
function to derive the categorization variables.

``` r
avalcax_lookup <- exprs(
  ~PARAMCD,            ~condition,             ~AVALCAT1,
  "TDURD",             AVAL >= 90,          ">= 90 days",
  "TDURD", AVAL >= 30 & AVAL < 90, ">= 30 and < 90 days",
  "TDURD",              AVAL < 30,           "< 30 days",
  "TDOSE",            AVAL < 1000,           "< 1000 mg",
  "TDOSE",           AVAL >= 1000,          ">= 1000 mg",
  "TPDOSE",           AVAL < 1000,           "< 1000 mg",
  "TPDOSE",          AVAL >= 1000,          ">= 1000 mg"
)

adex <- adex %>%
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  )
```

### Assign `ASEQ`

The [admiral](https://pharmaverse.github.io/admiral/) function
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_obs_number.md)
can be used to derive `ASEQ`. An example call is:

``` r
adex <- derive_var_obs_number(
  adex,
  new_var = ASEQ,
  by_vars = get_admiral_option("subject_keys"),
  order = exprs(PARCAT1, ASTDT, VISIT, VISITNUM, EXSEQ, PARAMN),
  check_type = "error"
)
```

### Add the `ADSL` variables

If needed, the other `ADSL` variables can now be added:

``` r
adex <- adex %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = get_admiral_option("subject_keys")
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
[pharmaverse E2E
example](https://pharmaverse.github.io/examples/adam/adsl).

## Example Script

| ADaM | Sourcing Command          |
|------|---------------------------|
| ADEX | `use_ad_template("ADEX")` |
