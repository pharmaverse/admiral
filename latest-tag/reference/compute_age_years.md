# Compute Age in Years

Converts a set of age values from the specified time unit to years.

## Usage

``` r
compute_age_years(age, age_unit)
```

## Arguments

- age:

  The ages to convert.

  A numeric vector is expected.

  Default value

  :   none

- age_unit:

  Age unit.

  Either a string containing the time unit of all ages in `age` or a
  character vector containing the time units of each age in `age` is
  expected. Note that permitted values are cases insensitive (e.g.
  `"YEARS"` is treated the same as `"years"` and `"Years"`).

  Permitted values

  :   `"years"`, `"months"`, `"weeks"`, `"days"`, `"hours"`,
      `"minutes"`, `"seconds"`, `NA_character_`.

  Default value

  :   none

## Value

The ages contained in `age` converted to years.

## Details

Returns a numeric vector of ages in years as doubles. Note that passing
`NA_character_` as a unit will result in an `NA` value for the outputted
age. Also note, underlying computations assume an equal number of days
in each year (365.25).

## See also

Date/Time Computation Functions that returns a vector:
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_dtf.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_duration.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dtm.md)

## Examples

``` r
compute_age_years(
  age = c(240, 360, 480),
  age_unit = "MONTHS"
)
#> [1] 20 30 40

compute_age_years(
  age = c(10, 520, 3650, 1000),
  age_unit = c("YEARS", "WEEKS", "DAYS", NA_character_)
)
#> [1] 10.000000  9.965777  9.993155        NA
```
