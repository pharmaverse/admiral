# Compute Duration

Compute duration between two dates, e.g., duration of an adverse event,
relative day, age, ...

## Usage

``` r
compute_duration(
  start_date,
  end_date,
  in_unit = "days",
  out_unit = "days",
  floor_in = TRUE,
  add_one = TRUE,
  trunc_out = FALSE,
  type = "duration"
)
```

## Arguments

- start_date:

  The start date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Refer to
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dt.md)
  to obtain a vector of imputed dates.

  Default value

  :   none

- end_date:

  The end date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Refer to
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dt.md)
  to obtain a vector of imputed dates.

  Default value

  :   none

- in_unit:

  Input unit

  See floor_in and add_one parameter for details.

  Permitted Values (case-insensitive):

  For years: `"year"`, `"years"`, `"yr"`, `"yrs"`, `"y"`

  For months: `"month"`, `"months"`, `"mo"`, `"mos"`

  For days: `"day"`, `"days"`, `"d"`

  For hours: `"hour"`, `"hours"`, `"hr"`, `"hrs"`, `"h"`

  For minutes: `"minute"`, `"minutes"`, `"min"`, `"mins"`

  For seconds: `"second"`, `"seconds"`, `"sec"`, `"secs"`, `"s"`

  Default value

  :   `"days"`

- out_unit:

  Output unit

  The duration is derived in the specified unit

  Permitted Values (case-insensitive):

  For years: `"year"`, `"years"`, `"yr"`, `"yrs"`, `"y"`

  For months: `"month"`, `"months"`, `"mo"`, `"mos"`

  For weeks: `"week"`, `"weeks"`, `"wk"`, `"wks"`, `"w"`

  For days: `"day"`, `"days"`, `"d"`

  For hours: `"hour"`, `"hours"`, `"hr"`, `"hrs"`, `"h"`

  For minutes: `"minute"`, `"minutes"`, `"min"`, `"mins"`

  For seconds: `"second"`, `"seconds"`, `"sec"`, `"secs"`, `"s"`

  Default value

  :   `"days"`

- floor_in:

  Round down input dates?

  The input dates are round down with respect to the input unit, e.g.,
  if the input unit is 'days', the time of the input dates is ignored.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- add_one:

  Add one input unit?

  If the duration is non-negative, one input unit is added. i.e., the
  duration can not be zero.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- trunc_out:

  Return integer part

  The fractional part of the duration (in output unit) is removed, i.e.,
  the integer part is returned.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- type:

  lubridate duration type.

  See below for details.

  Permitted values

  :   `"duration"`, `"interval"`

  Default value

  :   `"duration"`

## Value

The duration between the two date in the specified unit

## Details

The output is a numeric vector providing the duration as time from start
to end date in the specified unit. If the end date is before the start
date, the duration is negative.

## Duration Type

The [lubridate](https://lubridate.tidyverse.org/) package calculates two
types of spans between two dates: duration and interval. While these
calculations are largely the same, when the unit of the time period is
month or year the result can be slightly different.

The difference arises from the ambiguity in the length of `"1 month"` or
`"1 year"`. Months may have 31, 30, 28, or 29 days, and years are 365
days and 366 during leap years. Durations and intervals help solve the
ambiguity in these measures.

The **interval** between `2000-02-01` and `2000-03-01` is `1` (i.e. one
month). The **duration** between these two dates is `0.95`, which
accounts for the fact that the year 2000 is a leap year, February has 29
days, and the average month length is `30.4375`, i.e.
`29 / 30.4375 = 0.95`.

For additional details, review the [lubridate time span reference
page](https://lubridate.tidyverse.org/reference/timespan.html).

## See also

[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md)

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_age_years.md),
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_dtf.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dtm.md)

## Examples

``` r
library(lubridate)
#> 
#> Attaching package: ‘lubridate’
#> The following objects are masked from ‘package:base’:
#> 
#>     date, intersect, setdiff, union

# Derive duration in days (integer), i.e., relative day
compute_duration(
  start_date = ymd_hms("2020-12-06T15:00:00"),
  end_date = ymd_hms("2020-12-24T08:15:00")
)
#> [1] 19

# Derive duration in days (float)
compute_duration(
  start_date = ymd_hms("2020-12-06T15:00:00"),
  end_date = ymd_hms("2020-12-24T08:15:00"),
  floor_in = FALSE,
  add_one = FALSE
)
#> [1] 17.71875

# Derive age in years
compute_duration(
  start_date = ymd("1984-09-06"),
  end_date = ymd("2020-02-24"),
  trunc_out = TRUE,
  out_unit = "years",
  add_one = FALSE
)
#> [1] 35

# Derive duration in hours
compute_duration(
  start_date = ymd_hms("2020-12-06T9:00:00"),
  end_date = ymd_hms("2020-12-06T13:30:00"),
  out_unit = "hours",
  floor_in = FALSE,
  add_one = FALSE,
)
#> [1] 4.5
```
