# Compute Body Surface Area (BSA)

Computes BSA from height and weight making use of the specified
derivation method

## Usage

``` r
compute_bsa(height = height, weight = weight, method)
```

## Arguments

- height:

  HEIGHT value

  It is expected that HEIGHT is in cm.

  Permitted values

  :   numeric vector

  Default value

  :   `height`

- weight:

  WEIGHT value

  It is expected that WEIGHT is in kg.

  Permitted values

  :   numeric vector

  Default value

  :   `weight`

- method:

  Derivation method to use:

  Mosteller: sqrt(height \* weight / 3600)

  DuBois-DuBois: 0.007184 \* height ^ 0.725 \* weight ^ 0.425

  Haycock: 0.024265 \* height ^ 0.3964 \* weight ^ 0.5378

  Gehan-George: 0.0235 \* height ^ 0.42246 \* weight ^ 0.51456

  Boyd: 0.0003207 \* (height ^ 0.3) \* (1000 \* weight) ^ (0.7285 -
  (0.0188 \* log10(1000 \* weight)))

  Fujimoto: 0.008883 \* height ^ 0.663 \* weight ^ 0.444

  Takahira: 0.007241 \* height ^ 0.725 \* weight ^ 0.425

  Permitted values

  :   character value

  Default value

  :   none

## Value

The BSA (Body Surface Area) in m^2.

## Details

Usually this computation function can not be used with `%>%`.

## See also

[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bsa.md)

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bmi.md),
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
# Derive BSA by the Mosteller method
compute_bsa(
  height = 170,
  weight = 75,
  method = "Mosteller"
)
#> [1] 1.881932

# Derive BSA by the DuBois & DuBois method
compute_bsa(
  height = c(170, 185),
  weight = c(75, 90),
  method = "DuBois-DuBois"
)
#> [1] 1.863558 2.141011
```
