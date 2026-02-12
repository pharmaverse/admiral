# Get Imputation Targets

Determines the imputation targets for date (see
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_date.md)
and time (see
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_time.md))
components.

## Usage

``` r
get_imputation_targets(partial, date_imputation = NULL, time_imputation = NULL)
```

## Arguments

- partial:

  A list of partial date/time components.

  Default value

  :   none

- date_imputation:

  The value to impute the day/month when a datepart is missing.

  A character value is expected, either as a

  - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for
    the 15th of June,

  - or as a keyword: `"first"`, `"mid"`, `"last"` to impute to the
    first/mid/last day/month.

  Default value

  :   none

- time_imputation:

  The value to impute the time when a timepart is missing.

  A character value is expected, either as a

  - format with hour, min and sec specified as `"hh:mm:ss"`: e.g.
    `"00:00:00"` for the start of the day,

  - or as a keyword: `"first"`,`"last"` to impute to the start/end of a
    day.

  Default value

  :   none

## Value

A list of imputation targets for date and (if applicable) time
components.

## Examples

``` r
# Get imputation targets for a date with 'first' date imputation
partial_date <- list(year = "2020", month = "03", day = NA_character_)
target_first_date <- admiral:::get_imputation_targets(partial_date,
  date_imputation = "first",
  time_imputation = NULL
)
print(target_first_date)
#> $year
#> [1] "0000"
#> 
#> $month
#> [1] "01"
#> 
#> $day
#> [1] "01"
#> 

# Get imputation targets for a datetime with 'first' date and time imputation
partial_datetime <- list(
  year = "2020",
  month = "03",
  day = NA_character_,
  hour = "12",
  minute = NA_character_,
  second = NA_character_
)
target_first_datetime <- admiral:::get_imputation_targets(partial_datetime,
  date_imputation = "first",
  time_imputation = "first"
)
print(target_first_datetime)
#> $year
#> [1] "0000"
#> 
#> $month
#> [1] "01"
#> 
#> $day
#> [1] "01"
#> 
#> $hour
#> [1] "00"
#> 
#> $minute
#> [1] "00"
#> 
#> $second
#> [1] "00"
#> 

# Get imputation targets for a datetime with 'last' date and time imputation
target_last_datetime <- admiral:::get_imputation_targets(partial_datetime,
  date_imputation = "last",
  time_imputation = "last"
)
print(target_last_datetime)
#> $year
#> [1] "9999"
#> 
#> $month
#> [1] "12"
#> 
#> $day
#> [1] "28"
#> 
#> $hour
#> [1] "23"
#> 
#> $minute
#> [1] "59"
#> 
#> $second
#> [1] "59"
#> 

# Get imputation targets for a date with custom date imputation '06-15'
target_custom_date <- admiral:::get_imputation_targets(partial_date,
  date_imputation = "06-15",
  time_imputation = NULL
)
print(target_custom_date)
#> $year
#> [1] "xxxx"
#> 
#> $month
#> [1] "06"
#> 
#> $day
#> [1] "15"
#> 

# Get imputation targets for a datetime with custom time imputation '12:34:56'
target_custom_time <- admiral:::get_imputation_targets(partial_datetime,
  date_imputation = "first",
  time_imputation = "12:34:56"
)
print(target_custom_time)
#> $year
#> [1] "0000"
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
#> [1] "34"
#> 
#> $second
#> [1] "56"
#> 
```
