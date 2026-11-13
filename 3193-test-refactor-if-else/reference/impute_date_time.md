# Impute Missing Values

Imputes missing values in partial date/time components using target
values.

## Usage

``` r
impute_date_time(partial, target)
```

## Arguments

- partial:

  A list of partial date/time components.

  Default value

  :   none

- target:

  A list of target values for imputation.

  Default value

  :   none

## Value

A list of imputed date/time components.

## Examples

``` r
# Impute missing values for date components
partial_date <- list(year = "2020", month = NA_character_, day = NA_character_)
target_date <- list(year = "2020", month = "01", day = "01")
imputed_date <- admiral:::impute_date_time(partial_date, target_date)
print(imputed_date)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "01"
#> 
#> $day
#> [1] "01"
#> 

# Impute missing values for datetime components
partial_datetime <- list(
  year = "2020", month = NA_character_, day = NA_character_,
  hour = "12", minute = NA_character_, second = NA_character_
)
target_datetime <- list(
  year = "2020", month = "01", day = "01",
  hour = "12", minute = "00", second = "00"
)
imputed_datetime <- admiral:::impute_date_time(
  partial_datetime, target_datetime
)
print(imputed_datetime)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "01"
#> 
#> $day
#> [1] "01"
#> 
#> $hour
#> [1] "12"
#> 
#> $minute
#> [1] "00"
#> 
#> $second
#> [1] "00"
#> 

# Impute missing values when some components are already present
partial_mixed <- list(year = "2020", month = "06", day = NA_character_)
target_mixed <- list(year = "2020", month = "01", day = "01")
imputed_mixed <- admiral:::impute_date_time(partial_mixed, target_mixed)
print(imputed_mixed)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "06"
#> 
#> $day
#> [1] "01"
#> 
```
