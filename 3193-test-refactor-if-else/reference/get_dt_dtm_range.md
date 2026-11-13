# Get Range of Partial Date / Datetime

Internal helper function to convert a character vector of (possibly
partial) dates (`dtc`) into complete dates based on a specified
imputation rule (`date_imputation`).

## Usage

``` r
get_dt_dtm_range(
  dtc,
  lower_bounds = NULL,
  upper_bounds = NULL,
  create_datetime
)
```

## Arguments

- dtc:

  A character vector of dates in ISO 8601 format (e.g., `"2022-12-15"`,
  `"2022-12"`, `"2022"`). Partial dates are allowed.

  Default value

  :   none

- lower_bounds:

  Lower bounds restricting the range

  The specified bounds restrict the lower limit of the returned range if
  it is within the range of the possible dates.

  Permitted values

  :   a list of dates, e.g.
      `list(ymd_hms("2021-07-01T04:03:01"), ymd_hms("2022-05-12T13:57:23"))`

  Default value

  :   `NULL`

- upper_bounds:

  Upper bounds restricting the range

  The specified bounds restrict the upper limit of the returned range if
  it is within the range of the possible dates.

  Permitted values

  :   a list of dates, e.g.
      `list(ymd_hms("2021-07-01T04:03:01"), ymd_hms("2022-05-12T13:57:23"))`

  Default value

  :   `NULL`

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
library(lubridate)
# Get Range from Partial Dates
dtc_dates <- c("2020-02-29", "2021-03")
admiral:::get_dt_dtm_range(dtc_dates, create_datetime = FALSE)
#> $lower
#> [1] "2020-02-29" "2021-03-01"
#> 
#> $upper
#> [1] "2020-02-29" "2021-03-31"
#> 

# Get Range from Partial Datetime
dtc_datetimes <- c("2020-02-29T12:00", "2021-03T14:30")
admiral:::get_dt_dtm_range(dtc_datetimes, create_datetime = TRUE)
#> $lower
#> [1] "2020-02-29T12:00:00" "2021-03-01T00:00:00"
#> 
#> $upper
#> [1] "2020-02-29T12:00:59" "2021-03-31T23:59:59"
#> 

# Get Range with Bounds
dtc_dates <- c("2020-02-29", "2021-03")
admiral:::get_dt_dtm_range(
  dtc_dates,
  lower_bounds = list(c(ymd("2020-01-01"), ymd("2021-03-05"))),
  upper_bounds = list(c(ymd("2020-12-31"), ymd("2021-03-25"))),
  create_datetime = FALSE
)
#> $lower
#> [1] "2020-02-29" "2021-03-05"
#> 
#> $upper
#> [1] "2020-02-29" "2021-03-25"
#> 

# Edge case: Return empty character vector for empty input
admiral:::get_dt_dtm_range(character(0), create_datetime = TRUE)
#> character(0)
```
