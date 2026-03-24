# Derive the Time Imputation Flag

Derive the time imputation flag (`*TMF`) comparing a date character
vector (`--DTC`) with a Datetime vector (`*DTM`).

## Usage

``` r
compute_tmf(dtc, dtm, ignore_seconds_flag = TRUE)
```

## Arguments

- dtc:

  The date character vector (`--DTC`).

  A character date is expected in a format like `yyyy-mm-ddThh:mm:ss`
  (partial or complete).

  Default value

  :   none

- dtm:

  The Date vector to compare (`*DTM`).

  A datetime object is expected.

  Default value

  :   none

- ignore_seconds_flag:

  ADaM IG states that given SDTM (`--DTC`) variable, if only hours and
  minutes are ever collected, and seconds are imputed in (`*DTM`) as 00,
  then it is not necessary to set (`*TMF`) to `"S"`.

  By default it is assumed that no seconds are collected and `*TMF`
  shouldn't be set to `"S"`. A user can set this to `FALSE` if seconds
  are collected.

  The default value of `ignore_seconds_flag` is set to `TRUE` in admiral
  1.4.0 and later.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

## Value

The time imputation flag (`*TMF`) (character value of `"H"`, `"M"` ,
`"S"` or `NA`)

## Details

Usually this computation function can not be used with `%>%`.

## See also

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/compute_age_years.md),
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/main/reference/compute_dtf.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/main/reference/compute_duration.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md),
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/main/reference/convert_xxtpt_to_hours.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dtm.md)

## Examples

``` r
library(lubridate)

compute_tmf(dtc = "2019-07-18T15:25", dtm = ymd_hm("2019-07-18T15:25"))
#> [1] NA
compute_tmf(dtc = "2019-07-18T15", dtm = ymd_hm("2019-07-18T15:25"))
#> [1] "M"
compute_tmf(dtc = "2019-07-18", dtm = ymd("2019-07-18"))
#> [1] "H"
compute_tmf(dtc = "2022-05--T00:00", dtm = ymd_hm("2022-05-15T23:59"))
#> [1] "H"
compute_tmf(dtc = "2022-05--T23:00", dtm = ymd_hm("2022-05-15T23:59"))
#> [1] "M"
compute_tmf(
  dtc = "2022-05--T23:59:00",
  dtm = ymd_hms("2022-05-15T23:59:59"),
  ignore_seconds_flag = FALSE
)
#> [1] "S"
```
