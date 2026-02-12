# Derive the Date Imputation Flag

Derive the date imputation flag (`*DTF`) comparing a date character
vector (`--DTC`) with a Date vector (`*DT`).

## Usage

``` r
compute_dtf(dtc, dt)
```

## Arguments

- dtc:

  The date character vector (`--DTC`).

  A character date is expected in a format like `yyyy-mm-ddThh:mm:ss`
  (partial or complete).

  Default value

  :   none

- dt:

  The Date vector to compare.

  A date object is expected.

  Default value

  :   none

## Value

The date imputation flag (`*DTF`) (character value of `"D"`, `"M"` ,
`"Y"` or `NA`)

## Details

Usually this computation function can not be used with `%>%`.

## See also

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_age_years.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_duration.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/impute_dtc_dtm.md)

## Examples

``` r
compute_dtf(dtc = "2019-07", dt = as.Date("2019-07-18"))
#> [1] "D"
compute_dtf(dtc = "2019", dt = as.Date("2019-07-18"))
#> [1] "M"
compute_dtf(dtc = "--06-01T00:00", dt = as.Date("2022-06-01"))
#> [1] "Y"
compute_dtf(dtc = "2022-06--T00:00", dt = as.Date("2022-06-01"))
#> [1] "D"
compute_dtf(dtc = "2022---01T00:00", dt = as.Date("2022-06-01"))
#> [1] "M"
compute_dtf(dtc = "2022----T00:00", dt = as.Date("2022-06-01"))
#> [1] "M"
```
