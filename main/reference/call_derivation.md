# Call a Single Derivation Multiple Times

Call a single derivation multiple times with some parameters/arguments
being fixed across iterations and others varying.

## Usage

``` r
call_derivation(dataset = NULL, derivation, variable_params, ...)
```

## Arguments

- dataset:

  Input dataset

  Default value

  :   `NULL`

- derivation:

  The derivation function to call

  A function that performs a specific derivation is expected. A
  derivation adds variables or observations to a dataset. The first
  argument of a derivation must expect a dataset and the derivation must
  return a dataset. All expected arguments for the derivation function
  must be provided through the
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
  objects passed to the `variable_params` and `...` arguments.

  Default value

  :   none

- variable_params:

  A `list` of function arguments that are different across iterations.
  Each set of function arguments must be created using
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md).

  Default value

  :   none

- ...:

  Any number of *named* function arguments that stay the same across
  iterations. If a function argument is specified both inside
  `variable_params` and `...` then the value in `variable_params`
  overwrites the one in `...`.

  @details

  It is also possible to pass functions from outside the `{admiral}`
  package to `call_derivation()`, e.g. an extension package function, or
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
  The only requirement for a function being passed to `derivation` is
  that it must take a dataset as its first argument and return a
  dataset.

  Default value

  :   none

## Value

The input dataset with additional records/variables added depending on
which `derivation` has been used.

## See also

[`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)
`call_derivation()`

Higher Order Functions:
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md),
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)

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
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT, TRTEDT),
    by_vars = exprs(USUBJID)
  )

## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
## one can add multiple variables in one go
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
#> # A tibble: 16 × 9
#>    STUDYID DOMAIN USUBJID AESTDTC    AEENDTC    TRTSDT     TRTEDT     ASTDT     
#>    <chr>   <chr>  <chr>   <chr>      <chr>      <date>     <date>     <date>    
#>  1 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  2 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  3 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#>  4 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#>  5 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  6 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  7 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  8 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  9 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#> 10 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#> 11 PILOT01 AE     16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08
#> 12 PILOT01 AE     16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08
#> 13 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 14 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 15 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 16 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> # ℹ 1 more variable: AENDT <date>

## The above call using `call_derivation()` is equivalent to the following
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
#> # A tibble: 16 × 9
#>    STUDYID DOMAIN USUBJID AESTDTC    AEENDTC    TRTSDT     TRTEDT     ASTDT     
#>    <chr>   <chr>  <chr>   <chr>      <chr>      <date>     <date>     <date>    
#>  1 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  2 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  3 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#>  4 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#>  5 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  6 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  7 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  8 PILOT01 AE     06-1384 2012-09-15 2012-09-29 2012-09-15 2012-09-24 2012-09-15
#>  9 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#> 10 PILOT01 AE     06-1384 2012-09-23 2012-09-29 2012-09-15 2012-09-24 2012-09-23
#> 11 PILOT01 AE     16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08
#> 12 PILOT01 AE     16-1298 2013-06-08 2013-07-06 2013-04-08 2013-06-28 2013-06-08
#> 13 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 14 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 15 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> 16 PILOT01 AE     16-1298 2013-04-22 2013-07-06 2013-04-08 2013-06-28 2013-04-22
#> # ℹ 1 more variable: AENDT <date>
```
