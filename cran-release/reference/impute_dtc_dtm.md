# Impute Partial Date(-time) Portion of a `--DTC` Variable

Imputation partial date/time portion of a `--DTC` variable. based on
user input.

## Usage

``` r
impute_dtc_dtm(
  dtc,
  highest_imputation = "h",
  date_imputation = "first",
  time_imputation = "first",
  min_dates = NULL,
  max_dates = NULL,
  preserve = FALSE
)
```

## Arguments

- dtc:

  The `--DTC` date to impute

  A character date is expected in a format like `yyyy-mm-dd` or
  `yyyy-mm-ddThh:mm:ss`. Trailing components can be omitted and `-` is a
  valid "missing" value for any component.

  Permitted values

  :   a character date variable

  Default value

  :   none

- highest_imputation:

  Highest imputation level

  The `highest_imputation` argument controls which components of the
  `--DTC` value are imputed if they are missing. All components up to
  the specified level are imputed.

  If a component at a higher level than the highest imputation level is
  missing, `NA_character_` is returned. For example, for
  `highest_imputation = "D"` `"2020"` results in `NA_character_` because
  the month is missing.

  If `"n"` is specified, no imputation is performed, i.e., if any
  component is missing, `NA_character_` is returned.

  If `"Y"` is specified, `date_imputation` should be `"first"` or
  `"last"` and `min_dates` or `max_dates` should be specified
  respectively. Otherwise, `NA_character_` is returned if the year
  component is missing.

  Permitted values

  :   `"Y"` (year, highest level), `"M"` (month), `"D"` (day), `"h"`
      (hour), `"m"` (minute), `"s"` (second), `"n"` (none, lowest level)

  Default value

  :   `"h"`

- date_imputation:

  The value to impute the day/month when a datepart is missing.

  A character value is expected.

  - If `highest_imputation` is `"M"`, month and day can be specified as
    `"mm-dd"`: e.g. `"06-15"` for the 15th of June

  - When `highest_imputation` is `"M"` or `"D"`, the following keywords
    are available: `"first"`, `"mid"`, `"last"` to impute to the
    first/mid/last day/month. If `"mid"` is specified, missing
    components are imputed as the middle of the possible range:

    - If both month and day are missing, they are imputed as `"06-30"`
      (middle of the year).

    - If only day is missing, it is imputed as `"15"` (middle of the
      month).

  The year can not be specified; for imputing the year `"first"` or
  `"last"` together with `min_dates` or `max_dates` argument can be used
  (see examples).

  Permitted values

  :   `"first"`, `"mid"`, `"last"`, or user-defined

  Default value

  :   `"first"`

- time_imputation:

  The value to impute the time when a timepart is missing.

  A character value is expected, either as a

  - format with hour, min and sec specified as `"hh:mm:ss"`: e.g.
    `"00:00:00"` for the start of the day,

  - or as a keyword: `"first"`,`"last"` to impute to the start/end of a
    day.

  The argument is ignored if `highest_imputation = "n"`.

  Permitted values

  :   `"first"`, `"last"`, or user-defined

  Default value

  :   `"first"`

- min_dates:

  Minimum dates

  A list of dates is expected. It is ensured that the imputed date is
  not before any of the specified dates, e.g., that the imputed adverse
  event start date is not before the first treatment date. Only dates
  which are in the range of possible dates of the `dtc` value are
  considered. The possible dates are defined by the missing parts of the
  `dtc` date (see example below). This ensures that the non-missing
  parts of the `dtc` date are not changed. A date or date-time object is
  expected. For example

      impute_dtc_dtm(
        "2020-11",
        min_dates = list(
         ymd_hm("2020-12-06T12:12"),
         ymd_hm("2020-11-11T11:11")
        ),
        highest_imputation = "M"
      )

  returns `"2020-11-11T11:11:11"` because the possible dates for
  `"2020-11"` range from `"2020-11-01T00:00:00"` to
  `"2020-11-30T23:59:59"`. Therefore `"2020-12-06T12:12:12"` is ignored.
  Returning `"2020-12-06T12:12:12"` would have changed the month
  although it is not missing (in the `dtc` date).

  For date variables (not datetime) in the list the time is imputed to
  `"00:00:00"`. Specifying date variables makes sense only if the date
  is imputed. If only time is imputed, date variables do not affect the
  result.

  Permitted values

  :   a list of dates, e.g.
      `list(ymd_hms("2021-07-01T04:03:01"), ymd_hms("2022-05-12T13:57:23"))`

  Default value

  :   `NULL`

- max_dates:

  Maximum dates

  A list of dates is expected. It is ensured that the imputed date is
  not after any of the specified dates, e.g., that the imputed date is
  not after the data cut off date. Only dates which are in the range of
  possible dates are considered. A date or date-time object is expected.

  For date variables (not datetime) in the list the time is imputed to
  `"23:59:59"`. Specifying date variables makes sense only if the date
  is imputed. If only time is imputed, date variables do not affect the
  result.

  Permitted values

  :   a list of dates, e.g.
      `list(ymd_hms("2021-07-01T04:03:01"), ymd_hms("2022-05-12T13:57:23"))`

  Default value

  :   `NULL`

- preserve:

  Preserve lower level date/time part when higher order part is missing,
  e.g. preserve day if month is missing or preserve minute when hour is
  missing.

  For example `"2019---07"` would return `"2019-06-07` if
  `preserve = TRUE` (and `date_imputation = "mid"`).

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

## Value

A character vector

## Details

Usually this computation function can not be used with `%>%`.

## See also

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_age_years.md),
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_dtf.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_duration.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dt.md)

## Examples

``` r
library(lubridate)

dates <- c(
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019",
  "2019---07",
  ""
)

# No date imputation (highest_imputation defaulted to "h")
# Missing time part imputed with 00:00:00 portion by default
impute_dtc_dtm(dtc = dates)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:00" "2019-07-18T15:00:00"
#> [4] "2019-07-18T00:00:00" NA                    NA                   
#> [7] NA                    NA                    NA                   

# No date imputation (highest_imputation defaulted to "h")
# Missing time part imputed with 23:59:59 portion
impute_dtc_dtm(
  dtc = dates,
  time_imputation = "23:59:59"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:59" "2019-07-18T15:59:59"
#> [4] "2019-07-18T23:59:59" NA                    NA                   
#> [7] NA                    NA                    NA                   

# Same as above
impute_dtc_dtm(
  dtc = dates,
  time_imputation = "last"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:59" "2019-07-18T15:59:59"
#> [4] "2019-07-18T23:59:59" NA                    NA                   
#> [7] NA                    NA                    NA                   

# Impute to first day/month if date is partial
# Missing time part imputed with 00:00:00 portion by default
impute_dtc_dtm(
  dtc = dates,
  highest_imputation = "M"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:00" "2019-07-18T15:00:00"
#> [4] "2019-07-18T00:00:00" "2019-02-01T00:00:00" "2019-01-01T00:00:00"
#> [7] "2019-01-01T00:00:00" "2019-01-01T00:00:00" NA                   
# same as above
impute_dtc_dtm(
  dtc = dates,
  highest_imputation = "M",
  date_imputation = "01-01"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:00" "2019-07-18T15:00:00"
#> [4] "2019-07-18T00:00:00" "2019-02-01T00:00:00" "2019-01-01T00:00:00"
#> [7] "2019-01-01T00:00:00" "2019-01-01T00:00:00" NA                   

# Impute to last day/month if date is partial
# Missing time part imputed with 23:59:59 portion
impute_dtc_dtm(
  dtc = dates,
  date_imputation = "last",
  time_imputation = "last"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:59" "2019-07-18T15:59:59"
#> [4] "2019-07-18T23:59:59" NA                    NA                   
#> [7] NA                    NA                    NA                   

# Impute to mid day/month if date is partial
# Missing time part imputed with 00:00:00 portion by default
impute_dtc_dtm(
  dtc = dates,
  highest_imputation = "M",
  date_imputation = "mid"
)
#> [1] "2019-07-18T15:25:40" "2019-07-18T15:25:00" "2019-07-18T15:00:00"
#> [4] "2019-07-18T00:00:00" "2019-02-15T00:00:00" "2019-06-30T00:00:00"
#> [7] "2019-06-30T00:00:00" "2019-06-30T00:00:00" NA                   

# Impute a date and ensure that the imputed date is not before a list of
# minimum dates
impute_dtc_dtm(
  "2020-12",
  min_dates = list(
    ymd_hm("2020-12-06T12:12"),
    ymd_hm("2020-11-11T11:11")
  ),
  highest_imputation = "M"
)
#> [1] "2020-12-06T12:12:00"

# Impute completely missing dates (only possible if min_dates or max_dates is specified)
impute_dtc_dtm(
  c("2020-12", NA_character_),
  min_dates = list(
    ymd_hm("2020-12-06T12:12", "2020-01-01T01:01"),
    ymd_hm("2020-11-11T11:11", NA)
  ),
  highest_imputation = "Y"
)
#> [1] "2020-12-06T12:12:00" "2020-01-01T01:01:00"
```
