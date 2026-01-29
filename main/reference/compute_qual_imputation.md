# Function to Impute Values When Qualifier Exists in Character Result

Derive an imputed value

## Usage

``` r
compute_qual_imputation(character_value, imputation_type = 1, factor = 0)
```

## Arguments

- character_value:

  Character version of value to be imputed

  Default value

  :   none

- imputation_type:

  (default value=1) Valid Values: 1: Strip \<, \>, = and convert to
  numeric. 2: imputation_type=1 and if the character value contains a \<
  or \>, the number of of decimals associated with the character value
  is found and then a factor of 1/10^(number of decimals + 1) will be
  added/subtracted from the numeric value. If no decimals exists, a
  factor of 1/10 will be added/subtracted from the value.

  Default value

  :   `1`

- factor:

  Numeric value (default=0), when using `imputation_type` = 1, this
  value can be added or subtracted when the qualifier is removed.

  Default value

  :   `0`

## Value

The imputed value

## See also

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/main/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qtc.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/main/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/main/reference/transform_range.md)

## Examples

``` r
compute_qual_imputation("<40")
#> [1] 40
compute_qual_imputation(c("3", ">30.2"))
#> [1]  3.0 30.2
```
