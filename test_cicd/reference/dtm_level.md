# Create a `dtm_level` object

Create a `dtm_level` object

## Usage

``` r
dtm_level(level)
```

## Arguments

- level:

  Datetime level

  Permitted values

  :   `"Y"` (year, highest level), `"M"` (month), `"D"` (day), `"h"`
      (hour), `"m"` (minute), `"s"` (second, lowest level), `"n"` (none)

  Default value

  :   none

## Value

A `dtm_level` object

## Details

A `dtm_level` object is an ordered factor, i.e., two objects can be
compared.

## See also

Utilities used for date imputation:
[`dt_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dt_level.md),
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_date.md),
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_time.md),
[`get_partialdatetime()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_partialdatetime.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_imputed_dtc_dt.md),
[`restrict_imputed_dtc_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_imputed_dtc_dtm.md)

## Examples

``` r
# Create a dtm_level object with level "D" (day)
level_day <- admiral:::dtm_level("D")
print(level_day)
#> [1] D
#> Levels: n < s < m < h < D < M < Y

# Create a dtm_level object with level "h" (hour)
level_hour <- admiral:::dtm_level("h")
print(level_hour)
#> [1] h
#> Levels: n < s < m < h < D < M < Y

# Compare two dtm_level objects
level_day > level_hour # TRUE, because "D" is larger than "h".
#> [1] TRUE
```
