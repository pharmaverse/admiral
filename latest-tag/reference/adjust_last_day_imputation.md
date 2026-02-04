# Adjust Last Day Imputation

This functions adjusts the day of the imputed date to the last day the
month if the day was imputed. It should be called if
`date_imputation = "last"` was used for the date imputation as
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_imputation_target_date.md)
imputes the last day as `"28"`.

## Usage

``` r
adjust_last_day_imputation(imputed_dtc, partial)
```

## Arguments

- imputed_dtc:

  A character vector of imputed date/datetime strings.

  Default value

  :   none

- partial:

  A list of partial date/time components.

  Default value

  :   none

## Value

A character vector of adjusted date/datetime strings.

## Details

If the day component in `partial` is missing, the day (in `imputed_dtc`)
is adjusted to the last day of the month.

## Examples

``` r
# Adjust last day imputation for a date with an incomplete day
imputed_date <- "2021-03-28"
partial_date <- list(year = "2021", month = "03", day = NA_character_)
admiral:::adjust_last_day_imputation(imputed_date, partial_date)
#> [1] "2021-03-31"

# Adjust last day imputation for a datetime with missing day
imputed_datetime <- "2021-03-28T00:00:00"
partial_datetime <- list(
  year = "2021", month = "03", day = NA_character_,
  hour = "00", minute = "00", second = "00"
)
admiral:::adjust_last_day_imputation(imputed_datetime, partial_datetime)
#> [1] "2021-03-31T00:00:00"

# Adjust last day imputation for a date with known day
partial_date_known_day <- list(year = "2021", month = "03", day = "15")
adjusted_date_known_day <- admiral:::adjust_last_day_imputation(
  imputed_date,
  partial_date_known_day
)
print(adjusted_date_known_day)
#> [1] "2021-03-28"

# Adjust last day imputation for a datetime with known day
partial_datetime_known_day <- list(
  year = "2021", month = "03", day = "15",
  hour = "00", minute = "00", second = "00"
)
adjusted_datetime_known_day <- admiral:::adjust_last_day_imputation(
  imputed_datetime,
  partial_datetime_known_day
)
print(adjusted_datetime_known_day)
#> [1] "2021-03-28T00:00:00"
```
