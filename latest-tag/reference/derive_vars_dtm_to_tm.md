# Derive Time Variables from Datetime Variables

This function creates time variable(s) as output from datetime
variable(s)

## Usage

``` r
derive_vars_dtm_to_tm(dataset, source_vars)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `source_vars` argument are expected to
  be in the dataset.

  Default value

  :   none

- source_vars:

  A list of datetime variables created using
  [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md)
  from which time is to be extracted

  Default value

  :   none

## Value

A data frame containing the input dataset with the corresponding time
(`--TM`) variable(s) of all datetime variables (`--DTM`) specified in
`source_vars` with the correct name.

## Details

The names of the newly added variables are automatically set by
replacing the `--DTM` suffix of the `source_vars` with `--TM`. The
`--TM` variables are created using the `{hms}` package.

## See also

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_trtdurd.md),
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dt.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_dy.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

adcm <- tribble(
  ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM,
  "PAT01", "2012-02-25 23:41:10", "2012-02-28 19:03:00", "2013-02-25 23:32:16",
  "PAT01", "", "2012-02-28 19:00:00", "",
  "PAT01", "2017-02-25 23:00:02", "2013-02-25 19:00:15", "2014-02-25 19:00:56",
  "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:25:00", "2017-03-25 23:00:00",
  "PAT01", "2017-02-25 16:05:17", "2017-02-25 14:20:00", "2018-04-29 14:06:45",
) %>%
  mutate(
    TRTSDTM = as_datetime(TRTSDTM),
    ASTDTM = as_datetime(ASTDTM),
    AENDTM = as_datetime(AENDTM)
  )

adcm %>%
  derive_vars_dtm_to_tm(exprs(TRTSDTM)) %>%
  select(USUBJID, starts_with("TRT"), everything())
#> # A tibble: 5 × 5
#>   USUBJID TRTSDTM             TRTSTM   ASTDTM              AENDTM             
#>   <chr>   <dttm>              <time>   <dttm>              <dttm>             
#> 1 PAT01   2012-02-25 23:41:10 23:41:10 2012-02-28 19:03:00 2013-02-25 23:32:16
#> 2 PAT01   NA                        NA 2012-02-28 19:00:00 NA                 
#> 3 PAT01   2017-02-25 23:00:02 23:00:02 2013-02-25 19:00:15 2014-02-25 19:00:56
#> 4 PAT01   2017-02-25 16:00:00 16:00:00 2017-02-25 14:25:00 2017-03-25 23:00:00
#> 5 PAT01   2017-02-25 16:05:17 16:05:17 2017-02-25 14:20:00 2018-04-29 14:06:45

adcm %>%
  derive_vars_dtm_to_tm(exprs(TRTSDTM, ASTDTM, AENDTM)) %>%
  select(USUBJID, starts_with("TRT"), starts_with("AS"), starts_with("AE"))
#> # A tibble: 5 × 7
#>   USUBJID TRTSDTM             TRTSTM   ASTDTM              ASTTM   
#>   <chr>   <dttm>              <time>   <dttm>              <time>  
#> 1 PAT01   2012-02-25 23:41:10 23:41:10 2012-02-28 19:03:00 19:03:00
#> 2 PAT01   NA                        NA 2012-02-28 19:00:00 19:00:00
#> 3 PAT01   2017-02-25 23:00:02 23:00:02 2013-02-25 19:00:15 19:00:15
#> 4 PAT01   2017-02-25 16:00:00 16:00:00 2017-02-25 14:25:00 14:25:00
#> 5 PAT01   2017-02-25 16:05:17 16:05:17 2017-02-25 14:20:00 14:20:00
#> # ℹ 2 more variables: AENDTM <dttm>, AENTM <time>
```
