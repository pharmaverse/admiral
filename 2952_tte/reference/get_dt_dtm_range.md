# Get Range of Partial Date / Datetime

Internal helper function to convert a character vector of (possibly
partial) dates (`dtc`) into complete dates based on a specified
imputation rule (`date_imputation`).

## Usage

``` r
get_dt_dtm_range(dtc, create_datetime)
```

## Arguments

- dtc:

  A character vector of dates in ISO 8601 format (e.g., `"2022-12-15"`,
  `"2022-12"`, `"2022"`). Partial dates are allowed.

  Default value

  :   none

- create_datetime:

  return the range in datetime format.

  Default value

  :   none

## Value

A list containing two vectors of fully imputed dates in `"YYYY-MM-DD"`
or `"YYYY-MM-DDThh:mm:ss"` format - the lower and upper limit of the
range.

## Details

The functions replaces missing components in `dtc` with the earliest
(lower bound) and latest (upper bound) possible value. Missing year is
replaced with `"0000"` for the lower bound and `"9999"` for the upper
bound.

## Examples

``` r
# Get Range from Partial Dates
dtc_dates <- c("2020-02-29", "2021-03")
imputed_dates_first <- admiral:::get_dt_dtm_range(dtc_dates, create_datetime = FALSE)
print(imputed_dates_first)
#> $lower
#> [1] "2020-02-29" "2021-03-01"
#> 
#> $upper
#> [1] "2020-02-29" "2021-03-31"
#> 


# Get Range from Partial Datetime
dtc_datetimes <- c("2020-02-29T12:00", "2021-03T14:30")
imputed_datetimes_first <- admiral:::get_dt_dtm_range(dtc_datetimes, create_datetime = TRUE)
print(imputed_datetimes_first)
#> $lower
#> [1] "2020-02-29T12:00:00" "2021-03-01T00:00:00"
#> 
#> $upper
#> [1] "2020-02-29T12:00:59" "2021-03-31T23:59:59"
#> 

# Edge case: Return empty character vector for empty input
imputed_empty <- admiral:::get_dt_dtm_range(character(0), create_datetime = TRUE)
print(imputed_empty)
#> character(0)
```
