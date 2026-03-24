# Compute RR Interval From Heart Rate

Computes RR interval from heart rate.

## Usage

``` r
compute_rr(hr)
```

## Arguments

- hr:

  Heart rate

  A numeric vector is expected. It is expected that heart rate is
  measured in beats/min.

  Default value

  :   none

## Value

RR interval in ms: \$\$\frac{60000}{HR}\$\$

## Details

Usually this computation function can not be used with `%>%`.

## See also

[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md)

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/main/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation_dec.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/main/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/main/reference/transform_range.md)

## Examples

``` r
compute_rr(hr = 70.14)
#> [1] 855.432
```
