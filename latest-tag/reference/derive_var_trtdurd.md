# Derive Total Treatment Duration (Days)

Derives total treatment duration (days) (`TRTDURD`).

**Note:** This is a wrapper function for the more generic
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md).

## Usage

``` r
derive_var_trtdurd(dataset, start_date = TRTSDT, end_date = TRTEDT)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `start_date` and `end_date` arguments
  are expected to be in the dataset.

  Default value

  :   none

- start_date:

  The start date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `TRTSDT`

- end_date:

  The end date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `TRTEDT`

## Value

The input dataset with `TRTDURD` added

## Details

The total treatment duration is derived as the number of days from start
to end date plus one.

## See also

[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md)

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dy.md)

## Examples

``` r
library(tibble)
library(lubridate)

data <- tribble(
  ~TRTSDT, ~TRTEDT,
  ymd("2020-01-01"), ymd("2020-02-24")
)

derive_var_trtdurd(data)
#> # A tibble: 1 × 3
#>   TRTSDT     TRTEDT     TRTDURD
#>   <date>     <date>       <dbl>
#> 1 2020-01-01 2020-02-24      55
```
