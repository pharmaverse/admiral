# Derive Analysis Age

Derives analysis age (`AAGE`) and analysis age unit (`AAGEU`).

**Note:** This is a wrapper function for the more generic
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md).

## Usage

``` r
derive_vars_aage(
  dataset,
  start_date = BRTHDT,
  end_date = RANDDT,
  age_unit = "YEARS",
  type = "interval"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `start_date` and `end_date` arguments
  are expected to be in the dataset.

  Default value

  :   none

- start_date:

  The start date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `BRTHDT`

- end_date:

  The end date

  A date or date-time object is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   `RANDDT`

- age_unit:

  Age unit

  The age is derived in the specified unit

  Permitted values

  :   The values are considered case-insensitive.

      For years: `"year"`, `"years"`, `"yr"`, `"yrs"`, `"y"`

      For months: `"month"`, `"months"`, `"mo"`, `"mos"`

      For weeks: `"week"`, `"weeks"`, `"wk"`, `"wks"`, `"w"`

      For days: `"day"`, `"days"`, `"d"`

      For hours: `"hour"`, `"hours"`, `"hr"`, `"hrs"`, `"h"`

      For minutes: `"minute"`, `"minutes"`, `"min"`, `"mins"`

      For seconds: `"second"`, `"seconds"`, `"sec"`, `"secs"`, `"s"`

  Default value

  :   `"YEARS"`

- type:

  lubridate duration type

  See below for details.

  Default: `"interval"`

  Permitted Values: `"duration"`, `"interval"`

  Default value

  :   `"interval"`

## Value

The input dataset with `AAGE` and `AAGEU` added

## Details

The duration is derived as time from start to end date in the specified
output unit. If the end date is before the start date, the duration is
negative. The start and end date variable must be present in the
specified input dataset.

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

ADSL Functions that returns variable appended to dataset:
[`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_age_years.md),
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_extreme_event.md),
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_period.md)

## Examples

``` r
library(tibble)
library(lubridate)

data <- tribble(
  ~BRTHDT, ~RANDDT,
  ymd("1984-09-06"), ymd("2020-02-24")
)

derive_vars_aage(data)
#> # A tibble: 1 × 4
#>   BRTHDT     RANDDT      AAGE AAGEU
#>   <date>     <date>     <dbl> <chr>
#> 1 1984-09-06 2020-02-24    35 YEARS
```
