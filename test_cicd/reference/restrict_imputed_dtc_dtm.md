# Restrict Imputed `--DTC` date to Minimum/Maximum Dates

Restrict Imputed `--DTC` date to Minimum/Maximum Dates

## Usage

``` r
restrict_imputed_dtc_dtm(dtc, imputed_dtc, min_dates, max_dates)
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

- imputed_dtc:

  The imputed `--DTC` date

  Default value

  :   none

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

## Value

- The last of the minimum dates (`min_dates`) which are in the range of
  the partial `--DTC` date (`dtc`)

- The first of the maximum dates (`max_dates`) which are in the range of
  the partial `--DTC` date (`dtc`)

- `imputed_dtc` if the partial `--DTC` date (`dtc`) is not in range of
  any of the minimum or maximum dates.

## See also

[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dtm.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dt.md)

Utilities used for date imputation:
[`dt_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dt_level.md),
[`dtm_level()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dtm_level.md),
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_date.md),
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_imputation_target_time.md),
[`get_partialdatetime()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_partialdatetime.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/restrict_imputed_dtc_dt.md)
