# Compute Mean Arterial Pressure (MAP)

Computes mean arterial pressure (MAP) based on diastolic and systolic
blood pressure. Optionally heart rate can be used as well.

## Usage

``` r
compute_map(diabp, sysbp, hr = NULL)
```

## Arguments

- diabp:

  Diastolic blood pressure

  A numeric vector is expected.

  Default value

  :   none

- sysbp:

  Systolic blood pressure

  A numeric vector is expected.

  Default value

  :   none

- hr:

  Heart rate

  A numeric vector or `NULL` is expected.

  Default value

  :   `NULL`

## Value

A numeric vector of MAP values

## Details

\$\$\frac{2DIABP + SYSBP}{3}\$\$ if it is based on diastolic and
systolic blood pressure and \$\$DIABP + 0.01 e^{4.14 - \frac{40.74}{HR}}
(SYSBP - DIABP)\$\$ if it is based on diastolic, systolic blood
pressure, and heart rate.

Usually this computation function can not be used with `%>%`.

## See also

[`derive_param_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_map.md)

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_framingham.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/transform_range.md)

## Examples

``` r
# Compute MAP based on diastolic and systolic blood pressure
compute_map(diabp = 51, sysbp = 121)
#> [1] 74.33333

# Compute MAP based on diastolic and systolic blood pressure and heart rate
compute_map(diabp = 51, sysbp = 121, hr = 59)
#> [1] 73.03907
```
