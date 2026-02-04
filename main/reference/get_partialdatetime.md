# Parse `--DTC` variable and Determine Components

Parse `--DTC` variable and Determine Components

## Usage

``` r
get_partialdatetime(dtc, create_datetime)
```

## Arguments

- dtc:

  The `'--DTC'` date to parse

  A character date is expected in a format like `yyyy-mm-dd` or
  `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
  valid value for any component.

  Default value

  :   none

- create_datetime:

  logical scalar. If `TRUE` returns Datetime components. If `FALSE`
  returns Date components.

  Default value

  :   none

## Value

A list of character vectors. The elements of the list are named "year",
"month", "day", "hour", "minute", and "second". Missing components are
set to `NA_character_`.

## Details

The function can be replaced by the parttime parser once it is
available.

## See also

[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dtm.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dt.md)

Utilities used for date imputation:
[`dt_level()`](https:/pharmaverse.github.io/admiral/main/reference/dt_level.md),
[`dtm_level()`](https:/pharmaverse.github.io/admiral/main/reference/dtm_level.md),
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/main/reference/get_imputation_target_date.md),
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/main/reference/get_imputation_target_time.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_imputed_dtc_dt.md),
[`restrict_imputed_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_imputed_dtc_dtm.md)

## Examples

``` r
# Datetime
# Get partial datetime components for a complete datetime string
dtc_complete <- admiral:::get_partialdatetime("2020-03-15T12:34:56", TRUE)
print(dtc_complete)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] "15"
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

# Get partial datetime components for a partial datetime string
dtc_partial <- admiral:::get_partialdatetime("2020-03-15T12:34", TRUE)
print(dtc_partial)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] "15"
#> 
#> $hour
#> [1] "12"
#> 
#> $minute
#> [1] "34"
#> 
#> $second
#> [1] NA
#> 

# Get partial datetime components for a date-only string
dtc_date_only <- admiral:::get_partialdatetime("2020-03-15", TRUE)
print(dtc_date_only)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] "15"
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

# Get partial datetime components for an incomplete year string
dtc_year_partial <- admiral:::get_partialdatetime("2020", TRUE)
print(dtc_year_partial)
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

# Date
# Get partial date components for a complete datetime string
dtc_complete <- admiral:::get_partialdatetime("2020-03-15T12:34:56", FALSE)
print(dtc_complete)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] "15"
#> 

# Get partial date components for a partial datetime string
dtc_partial <- admiral:::get_partialdatetime("2020-03-15T12:34", FALSE)
print(dtc_partial)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] "15"
#> 

# Get partial date components for a year and month only string
dtc_month_only <- admiral:::get_partialdatetime("2020-03", FALSE)
print(dtc_month_only)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] "03"
#> 
#> $day
#> [1] NA
#> 

# Get partial date components for an incomplete year string
dtc_year_partial <- admiral:::get_partialdatetime("2020", FALSE)
print(dtc_year_partial)
#> $year
#> [1] "2020"
#> 
#> $month
#> [1] NA
#> 
#> $day
#> [1] NA
#> 
```
