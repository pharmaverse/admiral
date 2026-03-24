# Derive Change from Baseline

Derive change from baseline (`CHG`) in a BDS dataset

## Usage

``` r
derive_var_chg(dataset)
```

## Arguments

- dataset:

  Input dataset `AVAL` and `BASE` are expected.

  Default value

  :   none

## Value

The input dataset with an additional column named `CHG`

## Details

Change from baseline is calculated by subtracting the baseline value
from the analysis value.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_base.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)

advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
  "P01",    "WEIGHT", 80,    "Y",    80,
  "P01",    "WEIGHT", 80.8,  NA,     80,
  "P01",    "WEIGHT", 81.4,  NA,     80,
  "P02",    "WEIGHT", 75.3,  "Y",    75.3,
  "P02",    "WEIGHT", 76,    NA,     75.3
)
derive_var_chg(advs)
#> # A tibble: 5 × 6
#>   USUBJID PARAMCD  AVAL ABLFL  BASE   CHG
#>   <chr>   <chr>   <dbl> <chr> <dbl> <dbl>
#> 1 P01     WEIGHT   80   Y      80   0    
#> 2 P01     WEIGHT   80.8 NA     80   0.800
#> 3 P01     WEIGHT   81.4 NA     80   1.40 
#> 4 P02     WEIGHT   75.3 Y      75.3 0    
#> 5 P02     WEIGHT   76   NA     75.3 0.700
```
