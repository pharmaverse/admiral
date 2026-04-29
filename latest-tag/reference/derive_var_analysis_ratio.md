# Derive Ratio Variable

Derives a ratio variable for a BDS dataset based on user specified
variables.

## Usage

``` r
derive_var_analysis_ratio(dataset, numer_var, denom_var, new_var = NULL)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `numer_var` and `denom_var` arguments
  are expected to be in the dataset.

  Default value

  :   none

- numer_var:

  Variable containing numeric values to be used in the numerator of the
  ratio calculation.

  Default value

  :   none

- denom_var:

  Variable containing numeric values to be used in the denominator of
  the ratio calculation.

  Default value

  :   none

- new_var:

  A user-defined variable that will be appended to the dataset. The
  default behavior will take the denominator variable and prefix it with
  `R2` and append to the dataset. Using this argument will override this
  default behavior.

  Default is `NULL`.

  Default value

  :   `NULL`

## Value

The input dataset with a ratio variable appended

## Details

A user wishing to calculate a Ratio to Baseline, `AVAL / BASE` will have
returned a new variable `R2BASE` that will be appended to the input
dataset. Ratio to Analysis Range Lower Limit `AVAL / ANRLO` will return
a new variable `R2ANRLO`, and Ratio to Analysis Range Upper Limit
`AVAL / ANRHI` will return a new variable `R2ANRLO`. Please note how the
denominator variable has the prefix `R2----`. A user can override the
default returned variables by using the `new_var` argument. Also, values
of 0 in the denominator will return `NA` in the derivation.

Note that `R2AyHI` and `R2AyLO` can also be derived using this function.

Reference CDISC ADaM Implementation Guide Version 1.1 Section 3.3.4
Analysis Parameter Variables for BDS Datasets

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_basetype_records.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)

data <- tribble(
  ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
  "P01", "ALT", 1, 27, 27, 6, 34,
  "P01", "ALT", 2, 41, 27, 6, 34,
  "P01", "ALT", 3, 17, 27, 6, 34,
  "P02", "ALB", 1, 38, 38, 33, 49,
  "P02", "ALB", 2, 39, 38, 33, 49,
  "P02", "ALB", 3, 37, 38, 33, 49
)

# Returns "R2" prefixed variables
data %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE) %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO) %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI)
#> # A tibble: 6 × 10
#>   USUBJID PARAMCD   SEQ  AVAL  BASE ANRLO ANRHI R2BASE R2ANRLO R2ANRHI
#>   <chr>   <chr>   <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl>
#> 1 P01     ALT         1    27    27     6    34  1        4.5    0.794
#> 2 P01     ALT         2    41    27     6    34  1.52     6.83   1.21 
#> 3 P01     ALT         3    17    27     6    34  0.630    2.83   0.5  
#> 4 P02     ALB         1    38    38    33    49  1        1.15   0.776
#> 5 P02     ALB         2    39    38    33    49  1.03     1.18   0.796
#> 6 P02     ALB         3    37    38    33    49  0.974    1.12   0.755

# Returns user-defined variables
data %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE, new_var = R01BASE) %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO, new_var = R01ANRLO) %>%
  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI, new_var = R01ANRHI)
#> # A tibble: 6 × 10
#>   USUBJID PARAMCD   SEQ  AVAL  BASE ANRLO ANRHI R01BASE R01ANRLO R01ANRHI
#>   <chr>   <chr>   <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>
#> 1 P01     ALT         1    27    27     6    34   1         4.5     0.794
#> 2 P01     ALT         2    41    27     6    34   1.52      6.83    1.21 
#> 3 P01     ALT         3    17    27     6    34   0.630     2.83    0.5  
#> 4 P02     ALB         1    38    38    33    49   1         1.15    0.776
#> 5 P02     ALB         2    39    38    33    49   1.03      1.18    0.796
#> 6 P02     ALB         3    37    38    33    49   0.974     1.12    0.755
```
