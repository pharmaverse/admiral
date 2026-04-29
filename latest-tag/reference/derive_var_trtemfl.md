# Derive Treatment-emergent Flag

Derive treatment emergent analysis flag (e.g., `TRTEMFL`).

## Usage

``` r
derive_var_trtemfl(
  dataset,
  new_var = TRTEMFL,
  start_date = ASTDTM,
  end_date = AENDTM,
  trt_start_date = TRTSDTM,
  trt_end_date = NULL,
  end_window = NULL,
  ignore_time_for_trt_end = TRUE,
  initial_intensity = NULL,
  intensity = NULL,
  group_var = NULL,
  subject_keys = get_admiral_option("subject_keys")
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by `start_date`, `end_date`, `trt_start_date`,
  `trt_end_date`, `initial_intensity`, and `intensity` are expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- new_var:

  New variable

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `TRTEMFL`

- start_date:

  Event start date

  Permitted values

  :   a date or datetime variable

  Default value

  :   `ASTDTM`

- end_date:

  Event end date

  Permitted values

  :   a date or datetime variable

  Default value

  :   `AENDTM`

- trt_start_date:

  Treatment start date

  Permitted values

  :   a date or datetime variable

  Default value

  :   `TRTSDTM`

- trt_end_date:

  Treatment end date

  Permitted values

  :   a date or datetime variable

  Default value

  :   `NULL`

- end_window:

  If the argument is specified (in 'days'), events starting more than
  the specified number of days after end of treatment, are not flagged.

  Permitted values

  :   a positive integer, e.g. `2` or `5`

  Default value

  :   `NULL`

- ignore_time_for_trt_end:

  If the argument is set to `TRUE`, the time part is ignored for
  checking if the event occurred more than `end_window` days after end
  of treatment.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- initial_intensity:

  Initial severity/intensity or toxicity

  `initial_intensity` is ignored when `group_var` is specified.

  If this argument is specified and `group_var` is `NULL`, events which
  start before treatment start and end after treatment start (or are
  ongoing) and worsened (i.e., the intensity is greater than the initial
  intensity), are flagged.

  The values of the specified variable must be comparable with the usual
  comparison operators. I.e., if the intensity is greater than the
  initial intensity `initial_intensity < intensity` must evaluate to
  `TRUE`.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- intensity:

  Severity/intensity or toxicity

  If the argument is specified, events which start before treatment
  start and end after treatment start (or are ongoing) and worsened
  (i.e., the intensity is greater than the initial intensity), are
  flagged.

  The values of the specified variable must be comparable with the usual
  comparison operators. I.e., if the intensity is greater than the
  initial intensity `initial_intensity < intensity` must evaluate to
  `TRUE`.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- group_var:

  Grouping variable

  If the argument is specified, it assumes that AEs are recorded as one
  episode of AE with multiple lines using a grouping variable.

  Events starting during treatment or before treatment and worsening
  afterward are flagged. Once an AE record in a group is flagged, all
  subsequent records in the treatment window are flagged regardless of
  severity.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- subject_keys:

  Variables to uniquely identify a subject.

  This argument is only used when `group_var` is specified.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

The input dataset with the variable specified by `new_var` added

## Details

For the derivation of the new variable the following cases are
considered in this order. The first case which applies, defines the
value of the variable.

- *not treated*: If `trt_start_date` is `NA`, it is set to
  `NA_character_`.

- *event before treatment*: If `end_date` is before `trt_start_date`
  (and `end_date` is not `NA`), it is set to `NA_character_`.

- *no event date*: If `start_date` is `NA`, it is set to `"Y"` as in
  such cases it is usually considered more conservative to assume the
  event was treatment-emergent.

- *event started during treatment*:

  - if `end_window` is not specified: if `start_date` is on or after
    `trt_start_date`, it is set to `"Y"`,

  - if `end_window` is specified: if `start_date` is on or after
    `trt_start_date` and `start_date` is on or before `trt_end_date` +
    `end_window` days, it is set to `"Y"`,

- *event started before treatment and (possibly) worsened on treatment*:

  - if `initial_intensity`, `intensity` is specified and `group_var` is
    not specified: if `initial_intensity < intensity` and `start_date`
    is before `trt_start_date` and `end_date` is on or after
    `trt_start_date` or `end_date` is `NA`, it is set to `"Y"`;

  - if `group_var` is specified: if `intensity` at treatment start \<
    `intensity` and `start_date` is after `trt_start_date` and
    `end_date` is on or after `trt_start_date` or `end_date` is `NA`, it
    is set to `"Y"`;

- Otherwise it is set to `NA_character_`.

The behavior of `derive_var_trtemfl()` is aligned with the proposed
treatment-emergent AE assignment in the following [PHUSE White
Paper](https://phuse.s3.eu-central-1.amazonaws.com/Deliverables/Safety+Analytics/WP-087+Recommended+Definition+of++Treatment-Emergent+Adverse+Events+in+Clinical+Trials+.pdf).
See the final example in the examples section below.

## See also

OCCDS Functions:
[`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_atc.md),
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)

## Examples

### Basic treatment-emergent flag

Derive `TRTEMFL` without considering treatment end and worsening

- For this basic example, all we are using are AE start/end dates and
  comparing those against treatment start date.

- If the AE started on or after treatment then we flag as
  treatment-emergent (e.g. records 5-7).

- If missing AE start date then we flag as treatment-emergent as worst
  case (e.g. records 8, 11 and 13), unless we know that the AE end date
  was before treatment so we can rule out this being treatment-emergent
  (e.g. record 12).

- Any not treated subject would not get their AEs flagged as
  treatment-emergent (e.g. records 14-16).

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)
    library(lubridate)

    adae <- tribble(
      ~USUBJID, ~ASTDT,            ~AENDT,            ~AEITOXGR, ~AETOXGR,
      # before treatment
      "1",      ymd("2021-12-13"), ymd("2021-12-15"), "1",       "1",
      "1",      ymd("2021-12-14"), ymd("2021-12-14"), "1",       "3",
      # starting before treatment and ending during treatment
      "1",      ymd("2021-12-30"), ymd("2022-01-14"), "1",       "3",
      "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",       "1",
      # starting during treatment
      "1",      ymd("2022-01-01"), ymd("2022-01-02"), "3",       "4",
      # after treatment
      "1",      ymd("2022-05-10"), ymd("2022-05-10"), "2",       "2",
      "1",      ymd("2022-05-11"), ymd("2022-05-11"), "2",       "2",
      # missing dates
      "1",      NA,                NA,                "3",       "4",
      "1",      ymd("2021-12-30"), NA,                "3",       "4",
      "1",      ymd("2021-12-31"), NA,                "3",       "3",
      "1",      NA,                ymd("2022-01-04"), "3",       "4",
      "1",      NA,                ymd("2021-12-24"), "3",       "4",
      "1",      NA,                ymd("2022-06-04"), "3",       "4",
      # without treatment
      "2",      NA,                ymd("2021-12-03"), "1",       "2",
      "2",      ymd("2021-12-01"), ymd("2021-12-03"), "1",       "2",
      "2",      ymd("2021-12-06"), NA,                "1",       "2"
    ) %>%
      mutate(
        STUDYID = "AB42",
        TRTSDT = if_else(USUBJID == "1", ymd("2022-01-01"), NA),
        TRTEDT = if_else(USUBJID == "1", ymd("2022-04-30"), NA)
      )

    derive_var_trtemfl(
      adae,
      start_date = ASTDT,
      end_date = AENDT,
      trt_start_date = TRTSDT
    ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, TRTEMFL)
    #> # A tibble: 16 × 5
    #>    USUBJID TRTSDT     ASTDT      AENDT      TRTEMFL
    #>    <chr>   <date>     <date>     <date>     <chr>
    #>  1 1       2022-01-01 2021-12-13 2021-12-15 <NA>
    #>  2 1       2022-01-01 2021-12-14 2021-12-14 <NA>
    #>  3 1       2022-01-01 2021-12-30 2022-01-14 <NA>
    #>  4 1       2022-01-01 2021-12-31 2022-01-01 <NA>
    #>  5 1       2022-01-01 2022-01-01 2022-01-02 Y
    #>  6 1       2022-01-01 2022-05-10 2022-05-10 Y
    #>  7 1       2022-01-01 2022-05-11 2022-05-11 Y
    #>  8 1       2022-01-01 NA         NA         Y
    #>  9 1       2022-01-01 2021-12-30 NA         <NA>
    #> 10 1       2022-01-01 2021-12-31 NA         <NA>
    #> 11 1       2022-01-01 NA         2022-01-04 Y
    #> 12 1       2022-01-01 NA         2021-12-24 <NA>
    #> 13 1       2022-01-01 NA         2022-06-04 Y
    #> 14 2       NA         NA         2021-12-03 <NA>
    #> 15 2       NA         2021-12-01 2021-12-03 <NA>
    #> 16 2       NA         2021-12-06 NA         <NA>   

### Considering treatment end date (`trt_end_date` and `end_window`)

Derive `TRTEMFL` taking a treatment end window into account

- In addition to the treatment-emergent checks explained in the above
  example, we now supply a treatment end date, `trt_end_date = TRTEDT`
  and an end window, `end_window = 10`. With these, any AE which started
  on or before treatment end date + 10 days is considered as
  treatment-emergent. Otherwise, those starting after the treatment end
  window are no longer flagged as treatment-emergent (e.g. record 7).

    derive_var_trtemfl(
      adae,
      start_date = ASTDT,
      end_date = AENDT,
      trt_start_date = TRTSDT,
      trt_end_date = TRTEDT,
      end_window = 10
    ) %>% select(USUBJID, TRTSDT, TRTEDT, ASTDT, AENDT, TRTEMFL)
    #> # A tibble: 16 × 6
    #>    USUBJID TRTSDT     TRTEDT     ASTDT      AENDT      TRTEMFL
    #>    <chr>   <date>     <date>     <date>     <date>     <chr>
    #>  1 1       2022-01-01 2022-04-30 2021-12-13 2021-12-15 <NA>
    #>  2 1       2022-01-01 2022-04-30 2021-12-14 2021-12-14 <NA>
    #>  3 1       2022-01-01 2022-04-30 2021-12-30 2022-01-14 <NA>
    #>  4 1       2022-01-01 2022-04-30 2021-12-31 2022-01-01 <NA>
    #>  5 1       2022-01-01 2022-04-30 2022-01-01 2022-01-02 Y
    #>  6 1       2022-01-01 2022-04-30 2022-05-10 2022-05-10 Y
    #>  7 1       2022-01-01 2022-04-30 2022-05-11 2022-05-11 <NA>
    #>  8 1       2022-01-01 2022-04-30 NA         NA         Y
    #>  9 1       2022-01-01 2022-04-30 2021-12-30 NA         <NA>
    #> 10 1       2022-01-01 2022-04-30 2021-12-31 NA         <NA>
    #> 11 1       2022-01-01 2022-04-30 NA         2022-01-04 Y
    #> 12 1       2022-01-01 2022-04-30 NA         2021-12-24 <NA>
    #> 13 1       2022-01-01 2022-04-30 NA         2022-06-04 Y
    #> 14 2       NA         NA         NA         2021-12-03 <NA>
    #> 15 2       NA         NA         2021-12-01 2021-12-03 <NA>
    #> 16 2       NA         NA         2021-12-06 NA         <NA>   

### Considering treatment worsening (`initial_intensity` and `intensity`)

Derive a new variable named `TRTEM2FL` taking worsening after treatment
start into account

- We also now start look at changes in intensity following treatment
  start using the `initial_intensity` and `intensity` arguments. This
  only impacts AEs starting before treatment, and ending on or after
  treatment (or with missing AE end date). We can additionally consider
  treatment-emergence for an AE that was ongoing at the start of
  treatment which may have worsened as a result of treatment, i.e. the
  most extreme intensity is greater than the initial intensity (e.g.
  records 3 and 9).

    derive_var_trtemfl(
      adae,
      new_var = TRTEM2FL,
      start_date = ASTDT,
      end_date = AENDT,
      trt_start_date = TRTSDT,
      trt_end_date = TRTEDT,
      end_window = 10,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR
    ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, AEITOXGR, AETOXGR, TRTEM2FL)
    #> # A tibble: 16 × 7
    #>    USUBJID TRTSDT     ASTDT      AENDT      AEITOXGR AETOXGR TRTEM2FL
    #>    <chr>   <date>     <date>     <date>     <chr>    <chr>   <chr>
    #>  1 1       2022-01-01 2021-12-13 2021-12-15 1        1       <NA>
    #>  2 1       2022-01-01 2021-12-14 2021-12-14 1        3       <NA>
    #>  3 1       2022-01-01 2021-12-30 2022-01-14 1        3       Y
    #>  4 1       2022-01-01 2021-12-31 2022-01-01 1        1       <NA>
    #>  5 1       2022-01-01 2022-01-01 2022-01-02 3        4       Y
    #>  6 1       2022-01-01 2022-05-10 2022-05-10 2        2       Y
    #>  7 1       2022-01-01 2022-05-11 2022-05-11 2        2       <NA>
    #>  8 1       2022-01-01 NA         NA         3        4       Y
    #>  9 1       2022-01-01 2021-12-30 NA         3        4       Y
    #> 10 1       2022-01-01 2021-12-31 NA         3        3       <NA>
    #> 11 1       2022-01-01 NA         2022-01-04 3        4       Y
    #> 12 1       2022-01-01 NA         2021-12-24 3        4       <NA>
    #> 13 1       2022-01-01 NA         2022-06-04 3        4       Y
    #> 14 2       NA         NA         2021-12-03 1        2       <NA>
    #> 15 2       NA         2021-12-01 2021-12-03 1        2       <NA>
    #> 16 2       NA         2021-12-06 NA         1        2       <NA>    

### Worsening when the same AE is collected over multiple records (`intensity` and `group_var`)

Derive `TRTEMFL` taking worsening after treatment into account within a
grouping variable

- Firstly, to understand which records correspond to the same AE, we
  need to supply a grouping variable (`group_var`). Then this example
  works in a similar way to the above one, but here we don't have an
  initial intensity so we have to use the intensity of the AE at the
  time of treatment start. If an ongoing AE increases intensity after
  treatment start (i.e. worsens), then from that point on the records
  are considered treatment-emergent, unless after the treatment end
  window (e.g. records 4, 6 and 7).

    adae2 <- tribble(
      ~USUBJID, ~ASTDT,            ~AENDT,            ~AETOXGR, ~AEGRPID,
      # ongoing AE where intensity drops after treatment start
      "1",      ymd("2021-12-31"), ymd("2022-01-01"), "3",      "1",
      "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "1",
      # ongoing AE where intensity increases after treatment start
      "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",      "2",
      "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "2",
      # ongoing AE where intensity increases after treatment start and then drops
      "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",      "3",
      "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "3",
      "1",      ymd("2022-01-12"), ymd("2022-01-15"), "1",      "3"
    ) %>%
      mutate(
        STUDYID = "AB42",
        TRTSDT = if_else(USUBJID == "1", ymd("2022-01-01"), NA),
        TRTEDT = if_else(USUBJID == "1", ymd("2022-04-30"), NA)
      )

    derive_var_trtemfl(
      adae2,
      start_date = ASTDT,
      end_date = AENDT,
      trt_start_date = TRTSDT,
      trt_end_date = TRTEDT,
      end_window = 10,
      intensity = AETOXGR,
      group_var = AEGRPID
    ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, AETOXGR, AEGRPID, TRTEMFL)
    #> # A tibble: 7 × 7
    #>   USUBJID TRTSDT     ASTDT      AENDT      AETOXGR AEGRPID TRTEMFL
    #>   <chr>   <date>     <date>     <date>     <chr>   <chr>   <chr>
    #> 1 1       2022-01-01 2021-12-31 2022-01-01 3       1       <NA>
    #> 2 1       2022-01-01 2022-01-02 2022-01-11 2       1       <NA>
    #> 3 1       2022-01-01 2021-12-31 2022-01-01 1       2       <NA>
    #> 4 1       2022-01-01 2022-01-02 2022-01-11 2       2       Y
    #> 5 1       2022-01-01 2021-12-31 2022-01-01 1       3       <NA>
    #> 6 1       2022-01-01 2022-01-02 2022-01-11 2       3       Y
    #> 7 1       2022-01-01 2022-01-12 2022-01-15 1       3       Y      

### Further Examples from PHUSE White Paper

Here we present more cases (some new, some similar to the examples
above) which are aligned one-to-one with the scenarios in the [PHUSE
White
Paper](https://phuse.s3.eu-central-1.amazonaws.com/Deliverables/Safety+Analytics/WP-087+Recommended+Definition+of++Treatment-Emergent+Adverse+Events+in+Clinical+Trials+.pdf)

    adae3 <- tribble(
      ~USUBJID, ~TRTSDTM, ~TRTEDTM, ~ASTDTM, ~AENDTM, ~AEITOXGR, ~AETOXGR,
      # Patient 1: Pre-treatment AE
      "1", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
      # Patient 2: On-treatment AE
      "2", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
      # Patient 3: Pre-treatment AE, then on-treatment AE at same intensity
      "3", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
      "3", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
      # Patient 4: Pre-treatment AE, then on-treatment AE at wors. intensity
      "4", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
      "4", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "3",
      # Patient 5: Pre-treatment AE, then on-treatment AE at impr. intensity
      "5", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
      "5", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "1",
      # Patient 6: AE starting pre-treatment, continuing on-treatment, then 2nd AE at same intensity
      "6", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
      "6", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
      # Patient 7: AE starting pre-treatment, continuing on-treatment, then 2nd AE at wors. intensity
      "7", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
      "7", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "3",
      # Patient 8: AE starting pre-treatment, continuing on-treatment, then 2nd AE at impr. intensity
      "8", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
      "8", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "1",
      # Patient 9: AE starting pre-treatment, continuing on-treatment, and no change in intensity
      "9", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
      # Patient 10: AE starting pre-treatment, continuing on-treatment, and wors. intensity
      "10", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "4",
      # Patient 11: AE starting pre-treatment, continuing on-treatment, and impr. intensity
      "11", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "1",
      # Patient 12: AE starting pre-treatment, worsening, then improving
      "12", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "3", "2",
      # Patient 13: AE starting pre-treatment, improving, then worsening
      "13", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "1", "2",
    ) %>%
      mutate(
        ASTDTM = ymd(ASTDTM),
        AENDTM = ymd(AENDTM),
        TRTSDTM = ymd(TRTSDTM),
        TRTEDTM = ymd(TRTEDTM),
      )

    derive_var_trtemfl(
      adae3,
      new_var = TRTEMFL,
      trt_end_date = TRTEDTM,
      end_window = 0,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR,
      subject_keys = exprs(USUBJID)
    ) %>%
      select(USUBJID, TRTSDTM, TRTEDTM, ASTDTM, AENDTM, AEITOXGR, AETOXGR, TRTEMFL)
    #> # A tibble: 19 × 8
    #>    USUBJID TRTSDTM    TRTEDTM    ASTDTM     AENDTM     AEITOXGR AETOXGR TRTEMFL
    #>    <chr>   <date>     <date>     <date>     <date>     <chr>    <chr>   <chr>
    #>  1 1       2021-01-01 2021-12-31 2020-12-20 2020-12-21 2        2       <NA>
    #>  2 2       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        2       Y
    #>  3 3       2021-01-01 2021-12-31 2020-12-20 2020-12-21 2        2       <NA>
    #>  4 3       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        2       Y
    #>  5 4       2021-01-01 2021-12-31 2020-12-20 2020-12-21 2        2       <NA>
    #>  6 4       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        3       Y
    #>  7 5       2021-01-01 2021-12-31 2020-12-20 2020-12-21 2        2       <NA>
    #>  8 5       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        1       Y
    #>  9 6       2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        2       <NA>
    #> 10 6       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        2       Y
    #> 11 7       2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        2       <NA>
    #> 12 7       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        3       Y
    #> 13 8       2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        2       <NA>
    #> 14 8       2021-01-01 2021-12-31 2021-12-20 2021-12-21 2        1       Y
    #> 15 9       2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        2       <NA>
    #> 16 10      2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        4       Y
    #> 17 11      2021-01-01 2021-12-31 2020-12-23 2021-01-21 2        1       <NA>
    #> 18 12      2021-01-01 2021-12-31 2020-12-23 2021-01-21 3        2       <NA>
    #> 19 13      2021-01-01 2021-12-31 2020-12-23 2021-01-21 1        2       Y      
