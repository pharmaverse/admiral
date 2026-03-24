# Create dataset of single doses

Derives dataset of single dose from aggregate dose information. This may
be necessary when e.g. calculating last dose before an adverse event in
`ADAE` or deriving a total dose parameter in `ADEX` when
`EXDOSFRQ != ONCE`.

## Usage

``` r
create_single_dose_dataset(
  dataset,
  dose_freq = EXDOSFRQ,
  start_date = ASTDT,
  start_datetime = NULL,
  end_date = AENDT,
  end_datetime = NULL,
  lookup_table = dose_freq_lookup,
  lookup_column = CDISC_VALUE,
  nominal_time = NULL,
  keep_source_vars = expr_c(get_admiral_option("subject_keys"), dose_freq, start_date,
    start_datetime, end_date, end_datetime)
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `dose_freq`, `start_date`, and
  `end_date` arguments are expected to be in the dataset.

  Default value

  :   none

- dose_freq:

  The dose frequency

  The aggregate dosing frequency used for multiple doses in a row.

  Permitted values

  :   defined by lookup table.

  Default value

  :   `EXDOSFRQ`

- start_date:

  The start date

  A date object is expected. This object cannot contain `NA` values.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `ASTDT`

- start_datetime:

  The start date-time

  A date-time object is expected. This object cannot contain `NA`
  values.

  Refer to
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  to impute and derive a date-time from a date character vector to a
  date object.

  If the input dataset contains frequencies which refer to `DOSE_WINDOW`
  equals `"HOUR"` or `"MINUTE"`, the parameter must be specified.

  Default value

  :   `NULL`

- end_date:

  The end date

  A date or date-time object is expected. This object cannot contain
  `NA` values.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `AENDT`

- end_datetime:

  The end date-time

  A date-time object is expected. This object cannot contain `NA`
  values.

  Refer to
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  to impute and derive a date-time from a date character vector to a
  date object.

  If the input dataset contains frequencies which refer to `DOSE_WINDOW`
  equals `"HOUR"` or `"MINUTE"`, the parameter must be specified.

  Default value

  :   `NULL`

- lookup_table:

  The dose frequency value lookup table

  The table used to look up `dose_freq` values and determine the
  appropriate multiplier to be used for row generation. If a lookup
  table other than the default is used, it must have columns
  `DOSE_WINDOW`, `DOSE_COUNT`, and `CONVERSION_FACTOR`. The default
  table `dose_freq_lookup` is described in detail
  [here](https:/pharmaverse.github.io/admiral/main/reference/dose_freq_lookup.md).

  Permitted Values for `DOSE_WINDOW`: `"MINUTE"`, `"HOUR"`, `"DAY"`,
  `"WEEK"`, `"MONTH"`, `"YEAR"`

  Default value

  :   `dose_freq_lookup`

- lookup_column:

  The dose frequency value column in the lookup table

  The column of `lookup_table`.

  Default value

  :   `CDISC_VALUE`

- nominal_time:

  The nominal relative time from first dose (`NFRLT`)

  Used for PK analysis, this will be in hours and should be 0 for the
  first dose. It can be derived as `(VISITDY - 1) * 24` for example.
  This will be expanded as the single dose dataset is created. For
  example an `EXDOFRQ` of `"QD"` will result in the nominal_time being
  incremented by 24 hours for each expanded record.

  The value can be NULL if not needed.

  Default value

  :   `NULL`

- keep_source_vars:

  List of variables to be retained from source dataset

  This parameter can be specified if additional information is required
  in the output dataset. For example `EXTRT` for studies with more than
  one drug.

  Default value

  :   `expr_c(get_admiral_option("subject_keys"), dose_freq, start_date, start_datetime, end_date, end_datetime)`

## Value

The input dataset with a single dose per row.

## Details

Each aggregate dose row is split into multiple rows which each represent
a single dose.The number of completed dose periods between `start_date`
or `start_datetime` and `end_date` or `end_datetime` is calculated with
`compute_duration` and multiplied by `DOSE_COUNT`. For `DOSE_WINDOW`
values of `"WEEK"`, `"MONTH"`, and `"YEAR"`, `CONVERSION_FACTOR` is used
to convert into days the time object to be added to `start_date`.

Observations with dose frequency `"ONCE"` are copied to the output
dataset unchanged.

## See also

Creating auxiliary datasets:
[`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/main/reference/consolidate_metadata.md),
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_period_dataset.md),
[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)

## Examples

``` r
# Example with default lookup

library(lubridate)
library(stringr)
library(tibble)
library(dplyr)

data <- tribble(
  ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
  "STUDY01", "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
  ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
  "STUDY01", "P01", "Q3D", ymd("2021-01-08"), ymd_hms("2021-01-08 12:00:00"),
  ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
  "STUDY01", "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
  ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
)

create_single_dose_dataset(data)
#> # A tibble: 9 × 5
#>   STUDYID USUBJID EXDOSFRQ ASTDT      AENDT     
#>   <chr>   <chr>   <chr>    <date>     <date>    
#> 1 STUDY01 P01     ONCE     2021-01-01 2021-01-01
#> 2 STUDY01 P01     ONCE     2021-01-03 2021-01-03
#> 3 STUDY01 P01     ONCE     2021-01-05 2021-01-05
#> 4 STUDY01 P01     ONCE     2021-01-07 2021-01-07
#> 5 STUDY01 P01     ONCE     2021-01-08 2021-01-08
#> 6 STUDY01 P01     ONCE     2021-01-11 2021-01-11
#> 7 STUDY01 P01     ONCE     2021-01-14 2021-01-14
#> 8 STUDY01 P01     ONCE     2021-01-15 2021-01-15
#> 9 STUDY01 P01     ONCE     2021-01-29 2021-01-29

# Example with custom lookup

custom_lookup <- tribble(
  ~Value,   ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
  "Q30MIN", (1 / 30),    "MINUTE",                      1,
  "Q90MIN", (1 / 90),    "MINUTE",                      1
)

data <- tribble(
  ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
  "STUDY01", "P01", "Q30MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
  ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
  "STUDY02", "P02", "Q90MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
  ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00")
)

create_single_dose_dataset(data,
  lookup_table = custom_lookup,
  lookup_column = Value,
  start_datetime = ASTDTM,
  end_datetime = AENDTM
)
#> # A tibble: 6 × 7
#>   STUDYID USUBJID EXDOSFRQ ASTDT      ASTDTM              AENDT     
#>   <chr>   <chr>   <chr>    <date>     <dttm>              <date>    
#> 1 STUDY01 P01     ONCE     2021-01-01 2021-01-01 06:00:00 2021-01-01
#> 2 STUDY01 P01     ONCE     2021-01-01 2021-01-01 06:30:00 2021-01-01
#> 3 STUDY01 P01     ONCE     2021-01-01 2021-01-01 07:00:00 2021-01-01
#> 4 STUDY02 P02     ONCE     2021-01-01 2021-01-01 06:00:00 2021-01-01
#> 5 STUDY02 P02     ONCE     2021-01-01 2021-01-01 07:30:00 2021-01-01
#> 6 STUDY02 P02     ONCE     2021-01-01 2021-01-01 09:00:00 2021-01-01
#> # ℹ 1 more variable: AENDTM <dttm>
# Example with nominal time

data <- tribble(
  ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~NFRLT, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
  "STUDY01", "P01", "BID", 0, ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
  ymd("2021-01-07"), ymd_hms("2021-01-07 20:00:00"),
  "STUDY01", "P01", "BID", 168, ymd("2021-01-08"), ymd_hms("2021-01-08 08:00:00"),
  ymd("2021-01-14"), ymd_hms("2021-01-14 20:00:00"),
  "STUDY01", "P01", "BID", 336, ymd("2021-01-15"), ymd_hms("2021-01-15 08:00:00"),
  ymd("2021-01-29"), ymd_hms("2021-01-29 20:00:00")
)

create_single_dose_dataset(data,
  dose_freq = EXDOSFRQ,
  start_date = ASTDT,
  start_datetime = ASTDTM,
  end_date = AENDT,
  end_datetime = AENDTM,
  lookup_table = dose_freq_lookup,
  lookup_column = CDISC_VALUE,
  nominal_time = NFRLT,
  keep_source_vars = exprs(
    USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, NFRLT
  )
)
#> # A tibble: 58 × 7
#>    USUBJID EXDOSFRQ ASTDT      ASTDTM              AENDT     
#>    <chr>   <chr>    <date>     <dttm>              <date>    
#>  1 P01     ONCE     2021-01-01 2021-01-01 08:00:00 2021-01-01
#>  2 P01     ONCE     2021-01-01 2021-01-01 20:00:00 2021-01-01
#>  3 P01     ONCE     2021-01-02 2021-01-02 08:00:00 2021-01-02
#>  4 P01     ONCE     2021-01-02 2021-01-02 20:00:00 2021-01-02
#>  5 P01     ONCE     2021-01-03 2021-01-03 08:00:00 2021-01-03
#>  6 P01     ONCE     2021-01-03 2021-01-03 20:00:00 2021-01-03
#>  7 P01     ONCE     2021-01-04 2021-01-04 08:00:00 2021-01-04
#>  8 P01     ONCE     2021-01-04 2021-01-04 20:00:00 2021-01-04
#>  9 P01     ONCE     2021-01-05 2021-01-05 08:00:00 2021-01-05
#> 10 P01     ONCE     2021-01-05 2021-01-05 20:00:00 2021-01-05
#> # ℹ 48 more rows
#> # ℹ 2 more variables: AENDTM <dttm>, NFRLT <dbl>

# Example - derive a single dose dataset with imputations

# For either single drug administration records, or multiple drug administration
# records covering a range of dates, fill-in of missing treatment end datetime
# `EXENDTC` by substitution with an acceptable alternate, for example date of
# death, date of datacut may be required. This example shows the
# maximum possible number of single dose records to be derived. The example
# requires the date of datacut `DCUTDT` to be specified correctly, or
# if not appropriate to use `DCUTDT` as missing treatment end data and missing
# treatment end datetime could set equal to treatment start date and treatment
# start datetime. ADSL variables `DTHDT` and `DCUTDT` are preferred for
# imputation use.
#
# All available trial treatments are included, allowing multiple different
# last dose variables to be created in for example `use_ad_template("ADAE")`
# if required.

adsl <- tribble(
  ~STUDYID, ~USUBJID, ~DTHDT,
  "01", "1211", ymd("2013-01-14"),
  "01", "1083", ymd("2013-08-02"),
  "01", "1445", ymd("2014-11-01"),
  "01", "1015", NA,
  "01", "1023", NA
)

ex <- tribble(
  ~STUDYID, ~USUBJID, ~EXSEQ, ~EXTRT, ~EXDOSE, ~EXDOSU, ~EXDOSFRQ, ~EXSTDTC, ~EXENDTC,
  "01", "1015", 1, "PLAC", 0, "mg", "QD", "2014-01-02", "2014-01-16",
  "01", "1015", 2, "PLAC", 0, "mg", "QD", "2014-06-17", "2014-06-18",
  "01", "1015", 3, "PLAC", 0, "mg", "QD", "2014-06-19", NA_character_,
  "01", "1023", 1, "PLAC", 0, "mg", "QD", "2012-08-05", "2012-08-27",
  "01", "1023", 2, "PLAC", 0, "mg", "QD", "2012-08-28", "2012-09-01",
  "01", "1211", 1, "XANO", 54, "mg", "QD", "2012-11-15", "2012-11-28",
  "01", "1211", 2, "XANO", 54, "mg", "QD", "2012-11-29", NA_character_,
  "01", "1445", 1, "PLAC", 0, "mg", "QD", "2014-05-11", "2014-05-25",
  "01", "1445", 2, "PLAC", 0, "mg", "QD", "2014-05-26", "2014-11-01",
  "01", "1083", 1, "PLAC", 0, "mg", "QD", "2013-07-22", "2013-08-01"
)

adsl_death <- adsl %>%
  mutate(
    DTHDTM = convert_date_to_dtm(DTHDT),
    # Remove `DCUT` setup line below if ADSL `DCUTDT` is populated.
    DCUTDT = convert_dtc_to_dt("2015-03-06"), # Example only, enter date.
    DCUTDTM = convert_date_to_dtm(DCUTDT)
  )

# Select valid dose records, non-missing `EXSTDTC` and `EXDOSE`.
ex_mod <- ex %>%
  filter(!is.na(EXSTDTC) & !is.na(EXDOSE)) %>%
  derive_vars_merged(adsl_death, by_vars = get_admiral_option("subject_keys")) %>%
  # Example, set up missing `EXDOSFRQ` as QD daily dosing regime.
  # Replace with study dosing regime per trial treatment.
  mutate(EXDOSFRQ = if_else(is.na(EXDOSFRQ), "QD", EXDOSFRQ)) %>%
  # Create EXxxDTM variables and replace missing `EXENDTM`.
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    date_imputation = "first",
    time_imputation = "first",
    flag_imputation = "none",
  ) %>%
  derive_vars_dtm_to_dt(exprs(EXSTDTM)) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    # Maximum imputed treatment end date must not be not greater than
    # date of death or after the datacut date.
    max_dates = exprs(DTHDTM, DCUTDTM),
    date_imputation = "last",
    time_imputation = "last",
    flag_imputation = "none",
    highest_imputation = "Y",
  ) %>%
  derive_vars_dtm_to_dt(exprs(EXENDTM)) %>%
  # Select only unique values.
  # Removes duplicated records before final step.
  distinct(
    STUDYID, USUBJID, EXTRT, EXDOSE, EXDOSFRQ, DCUTDT, DTHDT, EXSTDT,
    EXSTDTM, EXENDT, EXENDTM, EXSTDTC, EXENDTC
  )

create_single_dose_dataset(
  ex_mod,
  start_date = EXSTDT,
  start_datetime = EXSTDTM,
  end_date = EXENDT,
  end_datetime = EXENDTM,
  keep_source_vars = exprs(
    STUDYID, USUBJID, EXTRT, EXDOSE, EXDOSFRQ,
    DCUTDT, EXSTDT, EXSTDTM, EXENDT, EXENDTM, EXSTDTC, EXENDTC
  )
)
#> # A tibble: 553 × 12
#>    STUDYID USUBJID EXTRT EXDOSE EXDOSFRQ DCUTDT     EXSTDT    
#>    <chr>   <chr>   <chr>  <dbl> <chr>    <date>     <date>    
#>  1 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-02
#>  2 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-03
#>  3 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-04
#>  4 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-05
#>  5 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-06
#>  6 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-07
#>  7 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-08
#>  8 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-09
#>  9 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-10
#> 10 01      1015    PLAC       0 ONCE     2015-03-06 2014-01-11
#> # ℹ 543 more rows
#> # ℹ 5 more variables: EXSTDTM <dttm>, EXENDT <date>, EXENDTM <dttm>,
#> #   EXSTDTC <chr>, EXENDTC <chr>
```
