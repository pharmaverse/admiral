# Compute Factor for Value Imputations When Character Value Contains \< or \>

Function to compute factor for value imputation when character value
contains \< or \>. The factor is calculated using the number of
decimals. If there are no decimals, the factor is 1, otherwise the
factor = 1/10^decimal place. For example, the factor for 100 = 1, the
factor for 5.4 = 1/10^1, the factor for 5.44 = 1/10^2. This results in
no additional false precision added to the value. This is an
intermediate function.

## Usage

``` r
compute_qual_imputation_dec(character_value_decimal)
```

## Arguments

- character_value_decimal:

  Character value to determine decimal precision

  Default value

  :   none

## Value

Decimal precision value to add or subtract

## Details

Derive an imputed value

## See also

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/transform_range.md)

## Examples

``` r
compute_qual_imputation_dec("<40.1")
#> [1] 0.1
compute_qual_imputation_dec(c("0.35", "1"))
#> [1] 0.01 1.00
```
