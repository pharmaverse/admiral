# Compute Corrected QT

Computes corrected QT using Bazett's, Fridericia's or Sagie's formula.

## Usage

``` r
compute_qtc(qt, rr, method)
```

## Arguments

- qt:

  QT interval

  A numeric vector is expected. It is expected that QT is measured in ms
  or msec.

  Default value

  :   none

- rr:

  RR interval

  A numeric vector is expected. It is expected that RR is measured in ms
  or msec.

  Default value

  :   none

- method:

  Method used to QT correction

  Permitted values

  :   `"Bazett"`, `"Fridericia"`, `"Sagie"`

  Default value

  :   none

## Value

QT interval in ms

## Details

Depending on the chosen `method` one of the following formulae is used.

*Bazett*: \$\$\frac{QT}{\sqrt{\frac{RR}{1000}}}\$\$

*Fridericia*: \$\$\frac{QT}{\sqrt\[3\]{\frac{RR}{1000}}}\$\$

*Sagie*: \$\$1000\left(\frac{QT}{1000} + 0.154\left(1 -
\frac{RR}{1000}\right)\right)\$\$

Usually this computation function can not be used with `%>%`.

## See also

[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md)

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/main/reference/compute_map.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/main/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/main/reference/transform_range.md)

## Examples

``` r
compute_qtc(qt = 350, rr = 857, method = "Bazett")
#> [1] 378.0747

compute_qtc(qt = 350, rr = 857, method = "Fridericia")
#> [1] 368.4748

compute_qtc(qt = 350, rr = 857, method = "Sagie")
#> [1] 372.022
```
