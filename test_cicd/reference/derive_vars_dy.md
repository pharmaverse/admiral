# Derive Relative Day Variables

Adds relative day variables (`*DY`) to the dataset, e.g., `ASTDY` and
`AENDY`.

## Usage

``` r
derive_vars_dy(dataset, reference_date, source_vars)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `reference_date` and `source_vars`
  arguments are expected to be in the dataset.

  Default value

  :   none

- reference_date:

  A date or date-time column, e.g., date of first treatment or date-time
  of last exposure to treatment.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   none

- source_vars:

  A list of datetime or date variables created using
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  from which dates are to be extracted. This can either be a list of
  date(time) variables or named `*DY` variables and corresponding
  `*DT(M)` variables e.g. `exprs(TRTSDTM, ASTDTM, AENDT)` or
  `exprs(TRTSDT, ASTDTM, AENDT, DEATHDY = DTHDT)`. If the source
  variable does not end in `*DT(M)`, a name for the resulting `*DY`
  variable must be provided.

  Default value

  :   none

## Value

The input dataset with `*DY` corresponding to the `*DTM` or `*DT` source
variable(s) added

## Details

The relative day is derived as number of days from the reference date to
the end date. If it is nonnegative, one is added. I.e., the relative day
of the reference date is 1. Unless a name is explicitly specified, the
name of the resulting relative day variable is generated from the source
variable name by replacing DT (or DTM as appropriate) with DY.

## See also

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_trtdurd.md),
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md)

## Examples

``` r
library(tibble)
library(lubridate)
library(dplyr, warn.conflicts = FALSE)

datain <- tribble(
  ~TRTSDTM, ~ASTDTM, ~AENDT,
  "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20"
) %>%
  mutate(
    TRTSDTM = as_datetime(TRTSDTM),
    ASTDTM = as_datetime(ASTDTM),
    AENDT = ymd(AENDT)
  )

derive_vars_dy(
  datain,
  reference_date = TRTSDTM,
  source_vars = exprs(TRTSDTM, ASTDTM, AENDT)
)
#> # A tibble: 1 × 6
#>   TRTSDTM             ASTDTM              AENDT      TRTSDY ASTDY AENDY
#>   <dttm>              <dttm>              <date>      <dbl> <dbl> <dbl>
#> 1 2014-01-17 23:59:59 2014-01-18 13:09:09 2014-01-20      1     2     4

# specifying name of new variables
datain <- tribble(
  ~TRTSDT, ~DTHDT,
  "2014-01-17", "2014-02-01"
) %>%
  mutate(
    TRTSDT = ymd(TRTSDT),
    DTHDT = ymd(DTHDT)
  )

derive_vars_dy(
  datain,
  reference_date = TRTSDT,
  source_vars = exprs(TRTSDT, DEATHDY = DTHDT)
)
#> # A tibble: 1 × 4
#>   TRTSDT     DTHDT      TRTSDY DEATHDY
#>   <date>     <date>      <dbl>   <dbl>
#> 1 2014-01-17 2014-02-01      1      16
```
