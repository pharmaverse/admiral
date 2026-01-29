# Create a `dt_level` object

Create a `dt_level` object

## Usage

``` r
dt_level(level)
```

## Arguments

- level:

  Date level

  Permitted values

  :   `"Y"` (year, highest level), `"M"` (month), `"D"` (day), `"n"`
      (none, lowest level)

  Default value

  :   none

## Value

A `dt_level` object

## Details

A `dt_level` object is an ordered factor, i.e., two objects can be
compared.

## See also

Utilities used for date imputation:
[`dtm_level()`](https:/pharmaverse.github.io/admiral/main/reference/dtm_level.md),
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/main/reference/get_imputation_target_date.md),
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/main/reference/get_imputation_target_time.md),
[`get_partialdatetime()`](https:/pharmaverse.github.io/admiral/main/reference/get_partialdatetime.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_imputed_dtc_dt.md),
[`restrict_imputed_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_imputed_dtc_dtm.md)

## Examples

``` r
# Create a dt_level object with level "D" (day)
level_day <- admiral:::dt_level("D")
print(level_day)
#> [1] D
#> Levels: n < D < M < Y

# Create a dt_level object with level "Y" (year)
level_year <- admiral:::dt_level("Y")
print(level_year)
#> [1] Y
#> Levels: n < D < M < Y

# Compare two dt_level objects
level_day > level_year # TRUE, because "Y" is larger than "D".
#> [1] FALSE
```
