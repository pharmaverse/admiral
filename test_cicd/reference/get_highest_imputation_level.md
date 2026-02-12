# Get Highest Imputation Level

Returns the
[`dt_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dt_level.md)
or
[`dtm_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dtm_level.md)
representation of the `highest_imputation` character value. The level
object allows comparisons of levels.

## Usage

``` r
get_highest_imputation_level(highest_imputation, create_datetime)
```

## Arguments

- highest_imputation:

  A character indicating the highest imputation level.

  Default value

  :   none

- create_datetime:

  A logical indicating whether datetime factors levels are required.

  Default value

  :   none

## Value

A
[`dt_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dt_level.md)
or
[`dtm_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dtm_level.md)
object representing the highest imputation level.

## Examples

``` r
# Get highest imputation level for date
highest_level_date <- admiral:::get_highest_imputation_level(
  highest_imputation = "Y",
  create_datetime = FALSE
)
print(highest_level_date)
#> [1] Y
#> Levels: n < D < M < Y

# Get highest imputation level for datetime
highest_level_datetime <- admiral:::get_highest_imputation_level(
  highest_imputation = "Y",
  create_datetime = TRUE
)
print(highest_level_datetime)
#> [1] Y
#> Levels: n < s < m < h < D < M < Y

# Get highest imputation level for date with month level
highest_level_month_date <- admiral:::get_highest_imputation_level(
  highest_imputation = "M",
  create_datetime = FALSE
)
print(highest_level_month_date)
#> [1] M
#> Levels: n < D < M < Y

# Get highest imputation level for datetime with hour level
highest_level_hour_datetime <- admiral:::get_highest_imputation_level(
  highest_imputation = "h",
  create_datetime = TRUE
)
print(highest_level_hour_datetime)
#> [1] h
#> Levels: n < s < m < h < D < M < Y
```
