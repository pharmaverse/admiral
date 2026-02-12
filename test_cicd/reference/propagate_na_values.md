# Propagate NA Values for datetime values

Propagates `NA` values through date/time components.

## Usage

``` r
propagate_na_values(partial)
```

## Arguments

- partial:

  A list of partial date/time components.

  Default value

  :   none

## Value

A list of date/time components with propagated `NA` values.

## Details

This function ensures that if a higher-order component (e.g., month) is
`NA`, all lower-order components (e.g., day, hour, etc.) are also set to
`NA`.

## Examples

``` r
# Propagate NA values through datetime components
partial_datetime <- list(
  year = "2020", month = NA_character_, day = "01",
  hour = "12", minute = NA_character_, second = "34"
)
propagated_datetime <- admiral:::propagate_na_values(partial_datetime)
print(propagated_datetime)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] NA
#> 
#> $day
#> [1] NA
#> 
#> $hour
#> [1] NA
#> 
#> $minute
#> [1] NA
#> 
#> $second
#> [1] NA
#> 

# Propagate NA values for datetime with missing higher order components
partial_missing <- list(
  year = NA_character_, month = "01", day = "01",
  hour = "12", minute = "00", second = "00"
)
propagated_missing <- admiral:::propagate_na_values(partial_missing)
print(propagated_missing)
#> $year
#> [1] NA
#> 
#> $month
#> [1] NA
#> 
#> $day
#> [1] NA
#> 
#> $hour
#> [1] NA
#> 
#> $minute
#> [1] NA
#> 
#> $second
#> [1] NA
#> 

partial_missing_date <- list(
  year = "2023", month = NA_character_, day = "01"
)
propagated_missing_date <- admiral:::propagate_na_values(partial_missing_date)
print(propagated_missing_date)
#> $year
#> [1] "2023"
#> 
#> $month
#> [1] NA
#> 
#> $day
#> [1] NA
#> 
```
