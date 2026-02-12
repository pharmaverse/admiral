# Derive Date Variables from Datetime Variables

This function creates date(s) as output from datetime variable(s)

## Usage

``` r
derive_vars_dtm_to_dt(dataset, source_vars)
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
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  from which dates are to be extracted

  Default value

  :   none

## Value

A data frame containing the input dataset with the corresponding date
(`--DT`) variable(s) of all datetime variables (`--DTM`) specified in
`source_vars.`

## See also

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_trtdurd.md),
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dy.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

adcm <- tribble(
  ~USUBJID, ~TRTSDTM,              ~ASTDTM,               ~AENDTM,
  "PAT01",  "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2012-02-25 23:00:00",
  "PAT01",  NA,                    "2012-02-28 19:00:00", NA,
  "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00",
  "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00",
  "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-04-29 14:00:00",
) %>%
  mutate(
    TRTSDTM = as_datetime(TRTSDTM),
    ASTDTM = as_datetime(ASTDTM),
    AENDTM = as_datetime(AENDTM)
  )

adcm %>%
  derive_vars_dtm_to_dt(exprs(TRTSDTM, ASTDTM, AENDTM)) %>%
  select(USUBJID, starts_with("TRT"), starts_with("AST"), starts_with("AEN"))
#> # A tibble: 5 × 7
#>   USUBJID TRTSDTM             TRTSDT     ASTDTM              ASTDT     
#>   <chr>   <dttm>              <date>     <dttm>              <date>    
#> 1 PAT01   2012-02-25 23:00:00 2012-02-25 2012-02-28 19:00:00 2012-02-28
#> 2 PAT01   NA                  NA         2012-02-28 19:00:00 2012-02-28
#> 3 PAT01   2017-02-25 23:00:00 2017-02-25 2013-02-25 19:00:00 2013-02-25
#> 4 PAT01   2017-02-25 16:00:00 2017-02-25 2017-02-25 14:00:00 2017-02-25
#> 5 PAT01   2017-02-25 16:00:00 2017-02-25 2017-02-25 14:00:00 2017-02-25
#> # ℹ 2 more variables: AENDTM <dttm>, AENDT <date>
```
