# Derive Duration

Derives duration between two dates, specified by the variables present
in input dataset e.g., duration of adverse events, relative day, age,
...

## Usage

``` r
derive_vars_duration(
  dataset,
  new_var,
  new_var_unit = NULL,
  start_date,
  end_date,
  in_unit = "days",
  out_unit = "DAYS",
  floor_in = TRUE,
  add_one = TRUE,
  trunc_out = FALSE,
  type = "duration"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `start_date` and `end_date` arguments
  are expected to be in the dataset.

  Default value

  :   none

- new_var:

  Name of variable to create

  Default value

  :   none

- new_var_unit:

  Name of the unit variable If the parameter is not specified, no
  variable for the unit is created.

  Default value

  :   `NULL`

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

The input dataset with the duration and unit variable added

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

[`compute_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_duration.md)

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_trtdurd.md),
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dy.md)

## Examples

``` r
library(lubridate)
library(tibble)

# Derive age in years
data <- tribble(
  ~USUBJID, ~BRTHDT, ~RANDDT,
  "P01", ymd("1984-09-06"), ymd("2020-02-24"),
  "P02", ymd("1985-01-01"), NA,
  "P03", NA, ymd("2021-03-10"),
  "P04", NA, NA
)

derive_vars_duration(data,
  new_var = AAGE,
  new_var_unit = AAGEU,
  start_date = BRTHDT,
  end_date = RANDDT,
  out_unit = "years",
  add_one = FALSE,
  trunc_out = TRUE
)
#> # A tibble: 4 × 5
#>   USUBJID BRTHDT     RANDDT      AAGE AAGEU
#>   <chr>   <date>     <date>     <dbl> <chr>
#> 1 P01     1984-09-06 2020-02-24    35 years
#> 2 P02     1985-01-01 NA            NA NA   
#> 3 P03     NA         2021-03-10    NA NA   
#> 4 P04     NA         NA            NA NA   

# Derive adverse event duration in days
data <- tribble(
  ~USUBJID, ~ASTDT, ~AENDT,
  "P01", ymd("2021-03-05"), ymd("2021-03-02"),
  "P02", ymd("2019-09-18"), ymd("2019-09-18"),
  "P03", ymd("1985-01-01"), NA,
  "P04", NA, NA
)

derive_vars_duration(data,
  new_var = ADURN,
  new_var_unit = ADURU,
  start_date = ASTDT,
  end_date = AENDT,
  out_unit = "days"
)
#> # A tibble: 4 × 5
#>   USUBJID ASTDT      AENDT      ADURN ADURU
#>   <chr>   <date>     <date>     <dbl> <chr>
#> 1 P01     2021-03-05 2021-03-02    -3 days 
#> 2 P02     2019-09-18 2019-09-18     1 days 
#> 3 P03     1985-01-01 NA            NA NA   
#> 4 P04     NA         NA            NA NA   

# Derive adverse event duration in minutes
data <- tribble(
  ~USUBJID, ~ADTM, ~TRTSDTM,
  "P01", ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-09T05:00:00"),
  "P02", ymd_hms("2019-11-11T10:30:00"), ymd_hms("2019-11-11T11:30:00"),
  "P03", ymd_hms("2019-11-11T00:00:00"), ymd_hms("2019-11-11T04:00:00"),
  "P04", NA, ymd_hms("2019-11-11T12:34:56"),
)

derive_vars_duration(data,
  new_var = ADURN,
  new_var_unit = ADURU,
  start_date = ADTM,
  end_date = TRTSDTM,
  in_unit = "minutes",
  out_unit = "minutes",
  add_one = FALSE
)
#> # A tibble: 4 × 5
#>   USUBJID ADTM                TRTSDTM             ADURN ADURU  
#>   <chr>   <dttm>              <dttm>              <dbl> <chr>  
#> 1 P01     2019-08-09 04:30:56 2019-08-09 05:00:00    30 minutes
#> 2 P02     2019-11-11 10:30:00 2019-11-11 11:30:00    60 minutes
#> 3 P03     2019-11-11 00:00:00 2019-11-11 04:00:00   240 minutes
#> 4 P04     NA                  2019-11-11 12:34:56    NA NA     

# Derive adverse event start time since last dose in hours
data <- tribble(
  ~USUBJID, ~ASTDTM, ~LDOSEDTM,
  "P01", ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-08T10:05:00"),
  "P02", ymd_hms("2019-11-11T23:59:59"), ymd_hms("2019-10-11T11:37:00"),
  "P03", ymd_hms("2019-11-11T00:00:00"), ymd_hms("2019-11-10T23:59:59"),
  "P04", ymd_hms("2019-11-11T12:34:56"), NA,
  "P05", NA, ymd_hms("2019-09-28T12:34:56")
)
derive_vars_duration(
  data,
  new_var = LDRELTM,
  new_var_unit = LDRELTMU,
  start_date = LDOSEDTM,
  end_date = ASTDTM,
  in_unit = "hours",
  out_unit = "hours",
  add_one = FALSE
)
#> # A tibble: 5 × 5
#>   USUBJID ASTDTM              LDOSEDTM            LDRELTM LDRELTMU
#>   <chr>   <dttm>              <dttm>                <dbl> <chr>   
#> 1 P01     2019-08-09 04:30:56 2019-08-08 10:05:00      18 hours   
#> 2 P02     2019-11-11 23:59:59 2019-10-11 11:37:00     756 hours   
#> 3 P03     2019-11-11 00:00:00 2019-11-10 23:59:59       1 hours   
#> 4 P04     2019-11-11 12:34:56 NA                       NA NA      
#> 5 P05     NA                  2019-09-28 12:34:56      NA NA      
```
