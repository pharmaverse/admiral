# Compute Body Mass Index (BMI)

Computes BMI from height and weight

## Usage

``` r
compute_bmi(height, weight)
```

## Arguments

- height:

  HEIGHT value

  It is expected that HEIGHT is in cm.

  Permitted values

  :   numeric vector

  Default value

  :   none

- weight:

  WEIGHT value

  It is expected that WEIGHT is in kg.

  Permitted values

  :   numeric vector

  Default value

  :   none

## Value

The BMI (Body Mass Index Area) in kg/m^2.

## Details

Usually this computation function can not be used with `%>%`.

## See also

[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bmi.md)

BDS-Findings Functions that returns a vector:
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/transform_range.md)

## Examples

``` r
compute_bmi(height = 170, weight = 75)
#> [1] 25.95156
```
