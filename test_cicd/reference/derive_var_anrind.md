# Derive Reference Range Indicator

Derive Reference Range Indicator

## Usage

``` r
derive_var_anrind(
  dataset,
  signif_dig = get_admiral_option("signif_digits"),
  use_a1hia1lo = FALSE
)
```

## Arguments

- dataset:

  Input dataset `ANRLO`, `ANRHI`, and `AVAL` are expected and if
  `use_a1hia1lo` is set to `TRUE`, `A1LO` and `A1H1` are expected as
  well.

  Default value

  :   none

- signif_dig:

  Number of significant digits to use when comparing values.

  Significant digits used to avoid floating point discrepancies when
  comparing numeric values. See blog: [How admiral handles floating
  points](https://pharmaverse.github.io/blog/posts/2023-10-30_floating_point/floating_point.html)

  Default value

  :   `get_admiral_option("signif_digits")`

- use_a1hia1lo:

  A logical value indicating whether to use `A1H1` and `A1LO` in the
  derivation of `ANRIND`.

  Default value

  :   `FALSE`

## Value

The input dataset with additional column `ANRIND`

## Details

In the case that `A1H1` and `A1LO` are to be used, `ANRIND` is set to:

- `"NORMAL"` if `AVAL` is greater or equal `ANRLO` and less than or
  equal `ANRHI`; or if `AVAL` is greater than or equal `ANRLO` and
  `ANRHI` is missing; or if `AVAL` is less than or equal `ANRHI` and
  `ANRLO` is missing

- `"LOW"` if `AVAL` is less than `ANRLO` and either `A1LO` is missing or
  `AVAL` is greater than or equal `A1LO`

- `"HIGH"` if `AVAL` is greater than `ANRHI` and either `A1HI` is
  missing or `AVAL` is less than or equal `A1HI`

- `"LOW LOW"` if `AVAL` is less than `A1LO`

- `"HIGH HIGH"` if `AVAL` is greater than `A1HI`

In the case that `A1H1` and `A1LO` are not to be used, `ANRIND` is set
to:

- `"NORMAL"` if `AVAL` is greater or equal `ANRLO` and less than or
  equal `ANRHI`; or if `AVAL` is greater than or equal `ANRLO` and
  `ANRHI` is missing; or if `AVAL` is less than or equal `ANRHI` and
  `ANRLO` is missing

- `"LOW"` if `AVAL` is less than `ANRLO`

- `"HIGH"` if `AVAL` is greater than `ANRHI`

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

vs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "P01",       "PUL",    70,     60,    100,    40,   110,
  "P01",       "PUL",    57,     60,    100,    40,   110,
  "P01",       "PUL",    60,     60,    100,    40,   110,
  "P01",     "DIABP",   102,     60,     80,    40,    90,
  "P02",       "PUL",   109,     60,    100,    40,   110,
  "P02",       "PUL",   100,     60,    100,    40,   110,
  "P02",     "DIABP",    80,     60,     80,    40,    90,
  "P03",       "PUL",    39,     60,    100,    40,   110,
  "P03",       "PUL",    40,     60,    100,    40,   110
)

vs %>% derive_var_anrind(use_a1hia1lo = TRUE)
#> # A tibble: 9 × 8
#>   USUBJID PARAMCD  AVAL ANRLO ANRHI  A1LO  A1HI ANRIND   
#>   <chr>   <chr>   <dbl> <dbl> <dbl> <dbl> <dbl> <chr>    
#> 1 P01     PUL        70    60   100    40   110 NORMAL   
#> 2 P01     PUL        57    60   100    40   110 LOW      
#> 3 P01     PUL        60    60   100    40   110 NORMAL   
#> 4 P01     DIABP     102    60    80    40    90 HIGH HIGH
#> 5 P02     PUL       109    60   100    40   110 HIGH     
#> 6 P02     PUL       100    60   100    40   110 NORMAL   
#> 7 P02     DIABP      80    60    80    40    90 NORMAL   
#> 8 P03     PUL        39    60   100    40   110 LOW LOW  
#> 9 P03     PUL        40    60   100    40   110 LOW      
vs %>% derive_var_anrind(use_a1hia1lo = FALSE)
#> # A tibble: 9 × 8
#>   USUBJID PARAMCD  AVAL ANRLO ANRHI  A1LO  A1HI ANRIND
#>   <chr>   <chr>   <dbl> <dbl> <dbl> <dbl> <dbl> <chr> 
#> 1 P01     PUL        70    60   100    40   110 NORMAL
#> 2 P01     PUL        57    60   100    40   110 LOW   
#> 3 P01     PUL        60    60   100    40   110 NORMAL
#> 4 P01     DIABP     102    60    80    40    90 HIGH  
#> 5 P02     PUL       109    60   100    40   110 HIGH  
#> 6 P02     PUL       100    60   100    40   110 NORMAL
#> 7 P02     DIABP      80    60    80    40    90 NORMAL
#> 8 P03     PUL        39    60   100    40   110 LOW   
#> 9 P03     PUL        40    60   100    40   110 LOW   
```
