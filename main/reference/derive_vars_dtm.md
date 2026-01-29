# Derive/Impute a Datetime from a Character Date

Derive a datetime object (`*DTM`) from a character date (`--DTC`). The
date and time can be imputed (see `date_imputation`/`time_imputation`
arguments) and the date/time imputation flag (`*DTF`, `*TMF`) can be
added.

## Usage

``` r
derive_vars_dtm(
  dataset,
  new_vars_prefix,
  dtc,
  highest_imputation = "h",
  date_imputation = "first",
  time_imputation = "first",
  flag_imputation = "auto",
  min_dates = NULL,
  max_dates = NULL,
  preserve = FALSE,
  ignore_seconds_flag = TRUE
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `dtc` argument are expected to be in
  the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- new_vars_prefix:

  Prefix used for the output variable(s).

  A character scalar is expected. For the date variable (`*DT`) is
  appended to the specified prefix, for the date imputation flag
  (`*DTF`), and for the time imputation flag (`*TMF`), i.e., for
  `new_vars_prefix = "AST"` the variables `ASTDT`, `ASTDTF`, and
  `ASTTMF` are created.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   none

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

- flag_imputation:

  Whether the date/time imputation flag(s) must also be derived.

  If `"both"` or `"date"` is specified, then date imputation flag is
  derived. If `"auto"` is specified and `highest_imputation` argument is
  greater than `"h"`, then date imputation flag is derived.

  If `"both"` or `"time"` is specified, then time imputation flag is
  derived. If `"auto"` is specified and `highest_imputation` argument is
  not `"n"`, then time imputation flag is derived.

  If `"none"` is specified, then no date or time imputation flag is
  derived.

  Please note that CDISC requirements dictate the need for a date/time
  imputation flag if any imputation is performed, so
  `flag_imputation = "none"` should only be used if the imputed variable
  is not part of the final ADaM dataset.

  Permitted values

  :   `"auto"`, `"date"`,`"time"`, `"both"` or `"none"`

  Default value

  :   `"auto"`

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

The input dataset with the datetime `*DTM` (and the date/time imputation
flag `*DTF`, `*TMF`) added.

## Details

In `{admiral}` we don't allow users to pick any single part of the
date/time to impute, we only enable to impute up to a highest level,
i.e. you couldn't choose to say impute months, but not days.

The presence of a `*DTF` variable is checked and the variable is not
derived if it already exists in the input dataset. However, if `*TMF`
already exists in the input dataset, a warning is issued and `*TMF` will
be overwritten.

## See also

[`vignette("imputation")`](https:/pharmaverse.github.io/admiral/main/articles/imputation.md)

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtdurd.md),
[`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)

## Examples

### Derive a datetime variable imputing time

In this example, we derive `ASTDTM` from `MHSTDTC`. Note that by default
the function imputes missing time components to `00` but doesn't impute
missing date components and automatically produces the time imputation
flag (`ASTTMF`).

    library(tibble)
    library(lubridate)

    mhdt <- tribble(
      ~MHSTDTC,
      "2019-07-18T15:25",
      "2019-07-18",
      "2019-02",
      "2019",
      "2019---07",
      ""
    )

    derive_vars_dtm(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC
    )
    #> # A tibble: 6 × 3
    #>   MHSTDTC            ASTDTM              ASTTMF
    #>   <chr>              <dttm>              <chr>
    #> 1 "2019-07-18T15:25" 2019-07-18 15:25:00 <NA>
    #> 2 "2019-07-18"       2019-07-18 00:00:00 H
    #> 3 "2019-02"          NA                  <NA>
    #> 4 "2019"             NA                  <NA>
    #> 5 "2019---07"        NA                  <NA>
    #> 6 ""                 NA                  <NA>  

### Impute to the latest (`date_imputation = "last"`)

In this example, we set `date_imputation = "last"` to get the last
month/day for partial dates. We also set `time_imputation = "last"`. The
function will use all or part of `23:59:59` for time imputation. Note
that `highest_imputation` must be at least `"D"` to perform date
imputation. Here we use `highest_imputation = "M"` to request imputation
of month and day (and time). Also note that two flag variables are
created. By default `ASTTMF` is set to `NA` if only seconds are imputed.
Set `ignore_seconds_flag = FALSE` to have the `"S"` flag for `ASTTMF`.

    derive_vars_dtm(
     mhdt,
     new_vars_prefix = "AST",
     dtc = MHSTDTC,
     date_imputation = "last",
     time_imputation = "last",
     highest_imputation = "M"
    )
    #> # A tibble: 6 × 4
    #>   MHSTDTC            ASTDTM              ASTDTF ASTTMF
    #>   <chr>              <dttm>              <chr>  <chr>
    #> 1 "2019-07-18T15:25" 2019-07-18 15:25:59 <NA>   <NA>
    #> 2 "2019-07-18"       2019-07-18 23:59:59 <NA>   H
    #> 3 "2019-02"          2019-02-28 23:59:59 D      H
    #> 4 "2019"             2019-12-31 23:59:59 M      H
    #> 5 "2019---07"        2019-12-31 23:59:59 M      H
    #> 6 ""                 NA                  <NA>   <NA>  

### Suppress imputation flags (`flag_imputation = "none"`)

In this example, we derive `ASTDTM` but suppress the `ASTTMF`. Note that
function appends missing `"hh:mm:ss"` to `ASTDTM`. The
`flag_imputation = "none"` call ensures no date/time imputation flag is
created. In practice, as per CDISC requirements this option can only be
selected if the imputed variable is not part of the final ADaM dataset.

    derive_vars_dtm(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      flag_imputation = "none"
    )
    #> # A tibble: 6 × 2
    #>   MHSTDTC            ASTDTM
    #>   <chr>              <dttm>
    #> 1 "2019-07-18T15:25" 2019-07-18 15:25:00
    #> 2 "2019-07-18"       2019-07-18 00:00:00
    #> 3 "2019-02"          NA
    #> 4 "2019"             NA
    #> 5 "2019---07"        NA
    #> 6 ""                 NA                 

### Avoid imputation after specified datetimes (`max_dates`)

In this example, we derive `AENDTM` where AE end date is imputed to the
last date. To ensure that the imputed date is not after the death or
data cut off date we can set `max_dates = exprs(DTHDT, DCUTDT)`. Note
two flag variables: `ASTDTF` and `ASTTMF` are created. Setting
`highest_imputation = "Y"` will allow for the missing `AEENDTC` record
to be imputed from `max_dates = exprs(DTHDT, DCUTDT)`.

    adae <- tribble(
       ~AEENDTC,             ~DTHDT,           ~DCUTDT,
       "2020-12", ymd("2020-12-26"), ymd("2020-12-24"),
       "2020-11", ymd("2020-12-06"), ymd("2020-12-24"),
              "", ymd("2020-12-06"), ymd("2020-12-24"),
    "2020-12-20", ymd("2020-12-06"), ymd("2020-12-24")
    )

    derive_vars_dtm(
      adae,
      dtc = AEENDTC,
      new_vars_prefix = "AEN",
      highest_imputation = "Y",
      date_imputation = "last",
      time_imputation = "last",
      max_dates = exprs(DTHDT, DCUTDT)
    )
    #> # A tibble: 4 × 6
    #>   AEENDTC      DTHDT      DCUTDT     AENDTM              AENDTF AENTMF
    #>   <chr>        <date>     <date>     <dttm>              <chr>  <chr>
    #> 1 "2020-12"    2020-12-26 2020-12-24 2020-12-24 23:59:59 D      H
    #> 2 "2020-11"    2020-12-06 2020-12-24 2020-11-30 23:59:59 D      H
    #> 3 ""           2020-12-06 2020-12-24 2020-12-06 23:59:59 Y      H
    #> 4 "2020-12-20" 2020-12-06 2020-12-24 2020-12-20 23:59:59 <NA>   H     

### Include `"S"` for time imputation flag (`ignore_seconds_flag`)

In this example, we set `ignore_seconds_flag = FALSE` to include `S` for
seconds in the `ASTTMF` variable. The default value of
`ignore_seconds_flag` is `TRUE` so the `"S"` is not normally displayed.
The ADaM IG states that given SDTM (`--DTC`) variable, if only hours and
minutes are ever collected, and seconds are imputed in (`*DTM`) as `00`,
then it is not necessary to set (`*TMF`) to `"S"`.

    mhdt <- tribble(
    ~MHSTDTC,
    "2019-07-18T15:25",
    "2019-07-18",
    "2019-02",
    "2019",
    "2019---07",
    ""
    )

    derive_vars_dtm(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      highest_imputation = "M",
      ignore_seconds_flag = FALSE
    )
    #> # A tibble: 6 × 4
    #>   MHSTDTC            ASTDTM              ASTDTF ASTTMF
    #>   <chr>              <dttm>              <chr>  <chr>
    #> 1 "2019-07-18T15:25" 2019-07-18 15:25:00 <NA>   S
    #> 2 "2019-07-18"       2019-07-18 00:00:00 <NA>   H
    #> 3 "2019-02"          2019-02-01 00:00:00 D      H
    #> 4 "2019"             2019-01-01 00:00:00 M      H
    #> 5 "2019---07"        2019-01-01 00:00:00 M      H
    #> 6 ""                 NA                  <NA>   <NA>  

### Preserve lower components if higher ones were imputed (`preserve`)

In this example, we impute dates as the middle month/day with
`date_imputation = "mid"` and impute time as last (`23:59:59`) with
`time_imputation = "last"`. We use the `preserve` argument to "preserve"
partial dates. For example, `"2019---18T15:-:05"`, will be displayed as
`"2019-06-18 15:59:05"` by setting `preserve = TRUE`.

    mhdt <- tribble(
    ~MHSTDTC,
    "2019-07-18T15:25",
    "2019---18T15:-:05",
    "2019-07-18",
    "2019-02",
    "2019",
    "2019---07",
    ""
    )

    derive_vars_dtm(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      highest_imputation = "M",
      date_imputation = "mid",
      time_imputation = "last",
      preserve = TRUE,
      ignore_seconds_flag = FALSE
    )
    #> # A tibble: 7 × 4
    #>   MHSTDTC             ASTDTM              ASTDTF ASTTMF
    #>   <chr>               <dttm>              <chr>  <chr>
    #> 1 "2019-07-18T15:25"  2019-07-18 15:25:59 <NA>   S
    #> 2 "2019---18T15:-:05" 2019-06-18 15:59:05 M      M
    #> 3 "2019-07-18"        2019-07-18 23:59:59 <NA>   H
    #> 4 "2019-02"           2019-02-15 23:59:59 D      H
    #> 5 "2019"              2019-06-30 23:59:59 M      H
    #> 6 "2019---07"         2019-06-07 23:59:59 M      H
    #> 7 ""                  NA                  <NA>   <NA>  

### Further examples

Further example usages of this function can be found in the
[`vignette("imputation")`](https:/pharmaverse.github.io/admiral/main/articles/imputation.md).
