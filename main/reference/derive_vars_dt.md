# Derive/Impute a Date from a Character Date

Derive a date (`*DT`) from a character date (`--DTC`). The date can be
imputed (see `date_imputation` argument) and the date imputation flag
(`*DTF`) can be added.

## Usage

``` r
derive_vars_dt(
  dataset,
  new_vars_prefix,
  dtc,
  highest_imputation = "n",
  date_imputation = "first",
  flag_imputation = "auto",
  min_dates = NULL,
  max_dates = NULL,
  preserve = FALSE
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
  appended to the specified prefix and for the date imputation flag
  (`*DTF`), i.e., for `new_vars_prefix = "AST"` the variables `ASTDT`
  and `ASTDTF` are created.

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

- flag_imputation:

  Whether the date imputation flag must also be derived.

  If `"auto"` is specified and `highest_imputation` argument is not
  `"n"`, then date imputation flag is derived.

  If `"date"` is specified, then date imputation flag is derived.

  If `"none"` is specified, then no date imputation flag is derived.

  Please note that CDISC requirements dictate the need for a date
  imputation flag if any imputation is performed, so
  `flag_imputation = "none"` should only be used if the imputed variable
  is not part of the final ADaM dataset.

  Permitted values

  :   `"auto"`, `"date"` or `"none"`

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

The input dataset with the date `*DT` (and the date imputation flag
`*DTF` if requested) added.

## Details

In `{admiral}` we don't allow users to pick any single part of the
date/time to impute, we only enable to impute up to a highest level,
i.e. you couldn't choose to say impute months, but not days.

The presence of a `*DTF` variable is checked and if it already exists in
the input dataset, a warning is issued and `*DTF` will be overwritten.

## See also

[`vignette("imputation")`](https:/pharmaverse.github.io/admiral/main/articles/imputation.md)

Date/Time Derivation Functions that returns variable appended to
dataset:
[`derive_var_trtdurd()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtdurd.md),
[`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
[`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_dt.md),
[`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_tm.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md),
[`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)

## Examples

### Derive a date variable without imputation

In this example, we derive `ASTDT` from `MHSTDTC` with no imputation
done for partial dates.

    library(tibble)
    library(lubridate)

    mhdt <- tribble(
      ~MHSTDTC,
      "2019-07-18T15:25:40",
      "2019-07-18T15:25",
      "2019-07-18",
      "2019-02",
      "2019",
      "2019---07",
      ""
    )

    derive_vars_dt(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC
    )
    #> # A tibble: 7 × 2
    #>   MHSTDTC               ASTDT
    #>   <chr>                 <date>
    #> 1 "2019-07-18T15:25:40" 2019-07-18
    #> 2 "2019-07-18T15:25"    2019-07-18
    #> 3 "2019-07-18"          2019-07-18
    #> 4 "2019-02"             NA
    #> 5 "2019"                NA
    #> 6 "2019---07"           NA
    #> 7 ""                    NA        

### Impute partial dates (`highest_imputation`)

Imputation is requested by the `highest_imputation` argument. Here
`highest_imputation = "M"` for month imputation is used, i.e. the
highest imputation done on a partial date is up to the month. By
default, missing date components are imputed to the first
day/month/year. A date imputation flag variable, `ASTDTF`, is
automatically created. The flag variable indicates if imputation was
done on the date.

    derive_vars_dt(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      highest_imputation = "M",
      date_imputation = "first"
    )
    #> # A tibble: 7 × 3
    #>   MHSTDTC               ASTDT      ASTDTF
    #>   <chr>                 <date>     <chr>
    #> 1 "2019-07-18T15:25:40" 2019-07-18 <NA>
    #> 2 "2019-07-18T15:25"    2019-07-18 <NA>
    #> 3 "2019-07-18"          2019-07-18 <NA>
    #> 4 "2019-02"             2019-02-01 D
    #> 5 "2019"                2019-01-01 M
    #> 6 "2019---07"           2019-01-01 M
    #> 7 ""                    NA         <NA>  

### Impute to the last day/month (`date_imputation = "last"`)

In this example, we derive `ADT` impute partial dates to last day/month,
i.e. `date_imputation = "last"`.

    qsdt <- tribble(
      ~QSDTC,
      "2019-07-18T15:25:40",
      "2019-07-18T15:25",
      "2019-07-18",
      "2019-02",
      "2019",
      "2019---07",
      ""
    )

    derive_vars_dt(
      qsdt,
      new_vars_prefix = "A",
      dtc = QSDTC,
      highest_imputation = "M",
      date_imputation = "last"
    )
    #> # A tibble: 7 × 3
    #>   QSDTC                 ADT        ADTF
    #>   <chr>                 <date>     <chr>
    #> 1 "2019-07-18T15:25:40" 2019-07-18 <NA>
    #> 2 "2019-07-18T15:25"    2019-07-18 <NA>
    #> 3 "2019-07-18"          2019-07-18 <NA>
    #> 4 "2019-02"             2019-02-28 D
    #> 5 "2019"                2019-12-31 M
    #> 6 "2019---07"           2019-12-31 M
    #> 7 ""                    NA         <NA> 

### Impute to the middle (`date_imputation = "mid"`) and suppress imputation flag (`flag_imputation = "none"`)

In this example, we will derive `TRTSDT` with date imputation flag
(`*DTF`) suppressed. Since `date_imputation = "mid"`, partial date
imputation will be set to June 30th for missing month and 15th for
missing day only. The `flag_imputation = "none"` call ensures no date
imputation flag is created. In practice, as per CDISC requirements this
option can only be selected if the imputed variable is not part of the
final ADaM dataset.

    exdt <- tribble(
      ~EXSTDTC,
      "2019-07-18T15:25:40",
      "2019-07-18T15:25",
      "2019-07-18",
      "2019-02",
      "2019",
      "2019---07",
      ""
    )
    derive_vars_dt(
      exdt,
      new_vars_prefix = "TRTS",
      dtc = EXSTDTC,
      highest_imputation = "M",
      date_imputation = "mid",
      flag_imputation = "none"
    )
    #> # A tibble: 7 × 2
    #>   EXSTDTC               TRTSDT
    #>   <chr>                 <date>
    #> 1 "2019-07-18T15:25:40" 2019-07-18
    #> 2 "2019-07-18T15:25"    2019-07-18
    #> 3 "2019-07-18"          2019-07-18
    #> 4 "2019-02"             2019-02-15
    #> 5 "2019"                2019-06-30
    #> 6 "2019---07"           2019-06-30
    #> 7 ""                    NA        

### Impute to a specific date (`date_imputation = "04-06"`)

In this example, we derive `ASTDT` with specific date imputation, i.e.
`date_imputation = "04-06"`. Note that day portion, `"-06"`, is used in
the imputation of the record with `"2019-02"`.

    derive_vars_dt(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      highest_imputation = "M",
      date_imputation = "04-06"
      )
    #> # A tibble: 7 × 3
    #>   MHSTDTC               ASTDT      ASTDTF
    #>   <chr>                 <date>     <chr>
    #> 1 "2019-07-18T15:25:40" 2019-07-18 <NA>
    #> 2 "2019-07-18T15:25"    2019-07-18 <NA>
    #> 3 "2019-07-18"          2019-07-18 <NA>
    #> 4 "2019-02"             2019-02-06 D
    #> 5 "2019"                2019-04-06 M
    #> 6 "2019---07"           2019-04-06 M
    #> 7 ""                    NA         <NA>  

### Applying a lower boundary to date imputation with (`min_dates`)

In this example, we derive `ASTDT` where `AESTDTC` is all partial dates
in need of imputation. Using `min_dates = exprs(TRTSDTM)`, we are
telling the function to apply the treatment start date `(TRTSDTM)` as a
lower boundary for imputation via the `min_dates` argument. This means:

- For partial dates that could potentially include `TRTSDTM` (case 1 &
  2), the imputed date is adjusted to `TRTSDTM`

- For partial dates that are entirely before `TRTSDTM` (case 3 & 4),
  standard imputation rules apply without adjustment

- For partial dates that are entirely after `TRTSDTM` (case 5), standard
  imputation rules apply

    adae <- tribble(
      ~case, ~AESTDTC, ~TRTSDTM,
      1, "2020-12", ymd_hms("2020-12-06T12:12:12"),
      2, "2020", ymd_hms("2020-12-06T12:12:12"),
      3, "2020-11", ymd_hms("2020-12-06T12:12:12"),
      4, "2020-01", ymd_hms("2020-12-06T12:12:12"),
      5, "2021-01", ymd_hms("2020-12-06T12:12:12")
    )

    derive_vars_dt(
      adae,
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "M",
      date_imputation = "first",
      min_dates = exprs(TRTSDTM)
      )
    #> # A tibble: 5 × 5
    #>    case AESTDTC TRTSDTM             ASTDT      ASTDTF
    #>   <dbl> <chr>   <dttm>              <date>     <chr>
    #> 1     1 2020-12 2020-12-06 12:12:12 2020-12-06 D
    #> 2     2 2020    2020-12-06 12:12:12 2020-12-06 M
    #> 3     3 2020-11 2020-12-06 12:12:12 2020-11-01 D
    #> 4     4 2020-01 2020-12-06 12:12:12 2020-01-01 D
    #> 5     5 2021-01 2020-12-06 12:12:12 2021-01-01 D     

### Applying an upper boundary to date imputation with (`max_dates`)

In this example, we derive `ASTDT` where `AESTDTC` is all partial dates
in need of imputation. Using `max_dates = exprs(TRTEDTM)`, we are
telling the function to apply the treatment end date `(TRTEDTM)` as an
upper boundary for imputation via the `max_dates` argument. This means:

- For partial dates that could potentially include `TRTEDTM` (case 1 &
  2), the imputed date is adjusted to `TRTEDTM`

- For partial dates that are entirely before `TRTEDTM` (case 3 & 4),
  standard imputation rules apply without adjustment

- For partial dates that are entirely after `TRTEDTM` (case 5), standard
  imputation rules apply

    adae <- tribble(
      ~case, ~AESTDTC, ~TRTSDTM, ~TRTEDTM,
      1, "2020-12", ymd_hms("2020-01-01T12:12:12"), ymd_hms("2020-12-20T23:59:59"),
      2, "2020", ymd_hms("2020-01-01T12:12:12"), ymd_hms("2020-12-20T23:59:59"),
      3, "2020-11", ymd_hms("2020-01-01T12:12:12"), ymd_hms("2020-12-20T23:59:59"),
      4, "2020-01", ymd_hms("2020-01-01T12:12:12"), ymd_hms("2020-12-20T23:59:59"),
      5, "2021-01", ymd_hms("2020-01-01T12:12:12"), ymd_hms("2020-12-20T23:59:59")
    )

    derive_vars_dt(
      adae,
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "M",
      date_imputation = "last",
      max_dates = exprs(TRTEDTM)
    )
    #> # A tibble: 5 × 6
    #>    case AESTDTC TRTSDTM             TRTEDTM             ASTDT      ASTDTF
    #>   <dbl> <chr>   <dttm>              <dttm>              <date>     <chr>
    #> 1     1 2020-12 2020-01-01 12:12:12 2020-12-20 23:59:59 2020-12-20 D
    #> 2     2 2020    2020-01-01 12:12:12 2020-12-20 23:59:59 2020-12-20 M
    #> 3     3 2020-11 2020-01-01 12:12:12 2020-12-20 23:59:59 2020-11-30 D
    #> 4     4 2020-01 2020-01-01 12:12:12 2020-12-20 23:59:59 2020-01-31 D
    #> 5     5 2021-01 2020-01-01 12:12:12 2020-12-20 23:59:59 2021-01-31 D     

### Preserve lower components if higher ones were imputed (`preserve`)

The `preserve` argument can be used to "preserve" information from the
partial dates. For example, `"2019---07"`, will be displayed as
`"2019-06-07"` rather than `"2019-06-30"` with `preserve = TRUE` and
`date_imputation = "mid"` .

    derive_vars_dt(
      mhdt,
      new_vars_prefix = "AST",
      dtc = MHSTDTC,
      highest_imputation = "M",
      date_imputation = "mid",
      preserve = TRUE
    )
    #> # A tibble: 7 × 3
    #>   MHSTDTC               ASTDT      ASTDTF
    #>   <chr>                 <date>     <chr>
    #> 1 "2019-07-18T15:25:40" 2019-07-18 <NA>
    #> 2 "2019-07-18T15:25"    2019-07-18 <NA>
    #> 3 "2019-07-18"          2019-07-18 <NA>
    #> 4 "2019-02"             2019-02-15 D
    #> 5 "2019"                2019-06-30 M
    #> 6 "2019---07"           2019-06-07 M
    #> 7 ""                    NA         <NA>  

### Further examples

Further example usages of this function can be found in the
[`vignette("imputation")`](https:/pharmaverse.github.io/admiral/main/articles/imputation.md).
