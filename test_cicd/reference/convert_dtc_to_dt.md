# Convert a Date Character Vector into a Date Object

Convert a date character vector (usually `--DTC`) into a Date vector
(usually `*DT`).

## Usage

``` r
convert_dtc_to_dt(
  dtc,
  highest_imputation = "n",
  date_imputation = "first",
  min_dates = NULL,
  max_dates = NULL,
  preserve = FALSE
)
```

## Arguments

- dtc:

  The –DTC date to convert.

  Permitted values

  :   a character date vector

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

  If `"n"` (none, lowest level) is specified no imputation is performed,
  i.e., if any component is missing, `NA_character_` is returned.

  If `"Y"` (year, highest level) is specified, `date_imputation` must be
  `"first"` or `"last"` and `min_dates` or `max_dates` must be specified
  respectively. Otherwise, an error is thrown.

  Permitted values

  :   `"Y"` (year, highest level), `"M"` (month), `"D"` (day), `"n"`
      (none, lowest level)

  Default value

  :   `"n"`

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
         ymd_hms("2020-12-06T12:12:12"),
         ymd_hms("2020-11-11T11:11:11")
        ),
        highest_imputation = "M"
      )

  returns `"2020-11-11T11:11:11"` because the possible dates for
  `"2020-11"` range from `"2020-11-01T00:00:00"` to
  `"2020-11-30T23:59:59"`. Therefore `"2020-12-06T12:12:12"` is ignored.
  Returning `"2020-12-06T12:12:12"` would have changed the month
  although it is not missing (in the `dtc` date).

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

  Permitted values

  :   a list of dates, e.g.
      `list(ymd_hms("2021-07-01T04:03:01"), ymd_hms("2022-05-12T13:57:23"))`

  Default value

  :   `NULL`

- preserve:

  Preserve day if month is missing and day is present

  For example `"2019---07"` would return `"2019-06-07` if
  `preserve = TRUE` (and `date_imputation = "MID"`).

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

## Value

a date object

## Details

Usually this computation function can not be used with `%>%`.

## See also

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_age_years.md),
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_dtf.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_duration.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dtm.md)

## Examples

``` r
convert_dtc_to_dt("2019-07-18")
#> [1] "2019-07-18"
convert_dtc_to_dt("2019-07")
#> [1] NA
```
