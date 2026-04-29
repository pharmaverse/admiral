# Create a Set of Parameters

Create a set of variable parameters/function arguments to be used in
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md).

## Usage

``` r
params(...)
```

## Arguments

- ...:

  One or more named arguments

  Default value

  :   none

## Value

An object of class `params`

## See also

[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md)

Other Advanced Functions:
[`list_tte_source_objects()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/list_tte_source_objects.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)

adsl <- tribble(
  ~STUDYID,   ~USUBJID,      ~TRTSDT,      ~TRTEDT,
  "PILOT01", "01-1307",           NA,           NA,
  "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
  "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
  "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
  "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
) %>%
  mutate(
    across(TRTSDT:TRTEDT, as.Date)
  )

ae <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
)

adae <- ae %>%
  select(USUBJID, AESTDTC, AEENDTC) %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT, TRTEDT),
    by_vars = exprs(USUBJID)
  )

## In order to derive both `ASTDT` and `AENDT` in `ADAE`, one can use `derive_vars_dt()`
adae %>%
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    date_imputation = "first",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  ) %>%
  derive_vars_dt(
    new_vars_prefix = "AEN",
    dtc = AEENDTC,
    date_imputation = "last",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  )
#> # A tibble: 16 × 7
#>    USUBJID AESTDTC    AEENDTC    TRTSDT     TRTEDT     ASTDT      AENDT     
#>    <chr>   <chr>      <chr>      <date>     <date>     <date>     <date>    
#>  1 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  2 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  3 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#>  4 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#>  5 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  6 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  7 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  8 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  9 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#> 10 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#> 11 16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08 2013-07-06
#> 12 16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08 2013-07-06
#> 13 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 14 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 15 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 16 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06


## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
## one can add multiple variables in one go.
## The function arguments which are different from a variable to another (e.g. `new_vars_prefix`,
## `dtc`, and `date_imputation`) are specified as a list of `params()` in the `variable_params`
## argument of `call_derivation()`. All other arguments which are common to all variables
## (e.g. `min_dates` and `max_dates`) are specified outside of `variable_params` (i.e. in `...`).
call_derivation(
  dataset = adae,
  derivation = derive_vars_dt,
  variable_params = list(
    params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
    params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
  ),
  min_dates = exprs(TRTSDT),
  max_dates = exprs(TRTEDT)
)
#> # A tibble: 16 × 7
#>    USUBJID AESTDTC    AEENDTC    TRTSDT     TRTEDT     ASTDT      AENDT     
#>    <chr>   <chr>      <chr>      <date>     <date>     <date>     <date>    
#>  1 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  2 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  3 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#>  4 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#>  5 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  6 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  7 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  8 06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15 2012-09-29
#>  9 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#> 10 06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23 2012-09-29
#> 11 16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08 2013-07-06
#> 12 16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08 2013-07-06
#> 13 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 14 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 15 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06
#> 16 16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22 2013-07-06

## The above call using `call_derivation()` is equivalent to the call using `derive_vars_dt()`
## to derive variables `ASTDT` and `AENDT` separately at the beginning.
```
