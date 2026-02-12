# Derive Nominal Relative Time from First Dose (NFRLT)

**\[experimental\]**

Derives nominal/planned time from first dose in hours by combining visit
day information with timepoint descriptions. The function converts
timepoint strings to hours using
[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md)
and adds them to the day-based offset. Optionally creates a
corresponding unit variable.

## Usage

``` r
derive_var_nfrlt(
  dataset,
  new_var = NFRLT,
  new_var_unit = NULL,
  out_unit = "HOURS",
  tpt_var = NULL,
  visit_day,
  first_dose_day = 1,
  treatment_duration = 0,
  range_method = "midpoint",
  set_values_to_na = NULL
)
```

## Arguments

- dataset:

  Input dataset containing visit day variable and optionally timepoint
  variable.

  Permitted values

  :   A data frame or tibble

  Default value

  :   none

- new_var:

  Name of the new variable to create (unquoted). Default is `NFRLT`.

  Permitted values

  :   Unquoted variable name

  Default value

  :   `NFRLT`

- new_var_unit:

  Name of the unit variable to create (unquoted). If specified, a
  character variable will be created containing the unit of time exactly
  as provided in `out_unit`. Common CDISC variables are `FRLTU` (First
  Dose Relative Time Unit) or `RRLTU` (Reference Relative Time Unit). If
  not specified, no unit variable is created.

  Permitted values

  :   Unquoted variable name (optional)

  Default value

  :   `NULL`

- out_unit:

  Unit of time for the output variable. Options are:

  - Days: "day", "days", "d"

  - Hours: "hour", "hours", "hr", "hrs", "h" (default: "hours")

  - Minutes: "minute", "minutes", "min", "mins"

  - Weeks: "week", "weeks", "wk", "wks", "w"

  Case-insensitive. The internal calculation is performed in hours, then
  converted to the specified unit. If `new_var_unit` is specified, it
  will contain the value exactly as provided by the user.

  Permitted values

  :   Character scalar (see options above)

  Default value

  :   `"HOURS"`

- tpt_var:

  Timepoint variable containing descriptions like "Pre-dose", "1H
  Post-dose", etc. (unquoted). If not provided or if the variable
  doesn't exist in the dataset, only the visit day offset is calculated
  (timepoint contribution is 0).

  Permitted values

  :   Unquoted variable name (optional)

  Default value

  :   `NULL`

- visit_day:

  Visit day variable (unquoted). This should be the planned/ nominal
  visit day (e.g., `VISITDY`). Records with `NA` in this variable will
  have NFRLT set to `NA`.

  Permitted values

  :   Unquoted variable name

  Default value

  :   none

- first_dose_day:

  The day number considered as the first dose day. Default is 1. For
  multiple-dose studies, this is typically Day 1.

  Permitted values

  :   Numeric scalar (positive integer)

  Default value

  :   `1`

- treatment_duration:

  Duration of treatment in hours. Can be either:

  - A numeric scalar (used for all records), or

  - An unquoted variable name from the dataset (e.g., `EXDUR`) where
    each record can have a different treatment duration

  Passed to
  [`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md).
  Must be non-negative. Default is 0 hours (for instantaneous treatments
  like oral medications).

  Permitted values

  :   Numeric scalar or unquoted variable name (non-negative)

  Default value

  :   `0`

- range_method:

  Method for converting time ranges to single values. Options are
  "midpoint" (default), "start", or "end". Passed to
  [`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md).
  For example, "0-6h" with midpoint returns 3, with start returns 0,
  with end returns 6.

  Permitted values

  :   Character scalar ("midpoint", "start", or "end")

  Default value

  :   `"midpoint"`

- set_values_to_na:

  An optional condition that marks derived NFRLT values as `NA`. For
  example, `set_values_to_na = VISIT == "UNSCHEDULED"` will set NFRLT to
  `NA` for all unscheduled visits. Can use any variables in the dataset.
  When `new_var_unit` is specified, the unit variable will also be set
  to `NA` for these records.

  Permitted values

  :   Condition (optional)

  Default value

  :   `NULL`

## Value

The input dataset with the new nominal relative time variable added, and
optionally the unit variable if `new_var_unit` is specified.

## Details

The nominal relative time is calculated as:

`NFRLT = (day_offset * 24 + timepoint_hours) * conversion_factor`

Where:

- `day_offset` is calculated from `visit_day` and `first_dose_day`,
  accounting for the absence of Day 0 in clinical trial convention

- `timepoint_hours` is derived from the timepoint description using
  [`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md),
  or 0 if `tpt_var` is not provided

- `conversion_factor` is:

  - 1 for "hours" (default)

  - 1/24 for "days"

  - 1/168 for "weeks" (1/24/7)

  - 60 for "minutes"

If `new_var_unit` is specified, a character variable is created
containing the value of `out_unit` exactly as provided by the user. For
example:

- `out_unit = "hours"` creates unit variable with value "hours"

- `out_unit = "HOURS"` creates unit variable with value "HOURS"

- `out_unit = "Days"` creates unit variable with value "Days"

- `NA` when the corresponding time value is `NA`

This matches the behavior of
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md)
and allows consistency when deriving multiple time variables.

**Handling "No Day 0":**

In clinical trials, day numbering typically follows the convention: ...,
Day -2, Day -1, Day 1, Day 2, ... (no Day 0). This function accounts for
this by adjusting the day offset when `visit_day` is negative and
`first_dose_day` is positive.

For example, with `first_dose_day = 1` and different output units:

- Day -1, `out_unit = "hours"` -\> -24 hours

- Day -1, `out_unit = "days"` -\> -1 day

- Day -1, `out_unit = "weeks"` -\> -0.1429 weeks

- Day -1, `out_unit = "minutes"` -\> -1440 minutes

- Day -7 -\> -168 hours, -7 days, -1 week, or -10080 minutes

- Day 1 -\> 0 (in any unit, first dose day)

- Day 8 -\> 168 hours, 7 days, 1 week, or 10080 minutes

With `first_dose_day = 7`:

- Day -1 -\> -168 hours, -7 days, -1 week, or -10080 minutes

- Day 1 -\> -144 hours, -6 days, -0.857 weeks, or -8640 minutes

- Day 6 -\> -24 hours, -1 day, -0.143 weeks, or -1440 minutes

- Day 7 -\> 0 (in any unit, first dose day)

**Common Use Cases:**

- **Single dose study**: Day 1 only, with samples at various timepoints
  (e.g., Pre-dose, 1H, 2H, 4H, 8H, 24H)

- **Multiple dose study**: Dosing on multiple days (e.g., Day 1, Day 8,
  Day 15) with samples around each dose

- **Screening visits**: Negative visit days (e.g., Day -14, Day -7)
  before first dose

- **Steady state study**: Multiple daily doses with sampling on specific
  days

- **Oral medications**: Use default `treatment_duration = 0` for
  instantaneous absorption

- **IV infusions**: Specify `treatment_duration` as infusion duration in
  hours (scalar) or as a variable name containing duration per record

- **Exposure records (EX)**: Can be called without `tpt_var` to derive
  NFRLT based only on visit day

- **Unscheduled visits**: Use `set_values_to_na` to set NFRLT to `NA`
  for unscheduled or early discontinuation visits

- **Variable treatment durations**: Use a variable name (e.g., `EXDUR`)
  when different subjects or visits have different treatment durations

- **Hours output**: Use `out_unit = "hours"` (default) for variables
  like `NFRLT` with `FRLTU`

- **Days output**: Use `out_unit = "days"` for variables like `NFRLTDY`
  with `FRLTU`

- **Weeks output**: Use `out_unit = "weeks"` for long-term studies with
  weekly dosing

- **Minutes output**: Use `out_unit = "minutes"` for very short-term PK
  studies or when minute precision is needed

- **CDISC compliance**: Use `new_var_unit = FRLTU` for first dose
  relative time or `new_var_unit = RRLTU` for reference relative time

- **Consistency with duration**: Use the same case for `out_unit` across
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md)
  and `derive_var_nfrlt()` to ensure unit variables match

**Important Notes:**

- The function assumes `visit_day` represents the nominal/planned day,
  not the actual study day

- Day numbering follows clinical trial convention with no Day 0

- For timepoints that span multiple days (e.g., "24H Post-dose"), ensure
  `visit_day` is set to the day when the sample was taken. For example,
  if dosing occurs on Day 3, a "24H Post-dose" sample taken on Day 4
  should have `visit_day = 4`.

- For crossover studies, consider deriving NFRLT separately per period

- `NA` values in `visit_day` will automatically result in `NA` for NFRLT
  (no need to use `set_values_to_na` for this case)

- `NA` values in `tpt_var` will result in `NA` for NFRLT

- `NA` values in the `treatment_duration` variable (if using a variable)
  will result in `NA` for NFRLT for those records

- Use `set_values_to_na` when you need to set NFRLT to `NA` based on
  other variables (e.g., `VISIT == "UNSCHEDULED"`), especially when
  `visit_day` is populated but should not be used for the NFRLT
  calculation

- If `tpt_var` is not provided or doesn't exist in the dataset,
  timepoint contribution is assumed to be 0 hours

- When using non-hour units, timepoint contributions are still
  calculated in hours first (e.g., "2H Post-dose" = 2 hours), then the
  entire result is converted to the specified unit

- The unit variable (if created) will contain the exact value provided
  in `out_unit`, preserving case and format

**Setting Special Values:**

If you need to set NFRLT to a specific value (e.g., 99999) for certain
visits instead of `NA`, use `set_values_to_na` first to set them to
`NA`, then use a subsequent
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) call to
replace those `NA` values:

    dataset %>%
      derive_var_nfrlt(
        ...,
        set_values_to_na = VISIT == "UNSCHEDULED"
      ) %>%
      mutate(NFRLT = if_else(is.na(NFRLT) & VISIT == "UNSCHEDULED", 99999, NFRLT))

## See also

[`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_xxtpt_to_hours.md),
[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_duration.md)

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_crit_flag.md)

## Examples

### Single dose study

Day 1 only with oral medication

    library(dplyr)
    library(tibble)

    adpc <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    1,        "1H Post-dose",
      "001",    1,        "2H Post-dose",
      "001",    1,        "4H Post-dose",
      "001",    1,        "24H Post-dose"
    )

    derive_var_nfrlt(
      adpc,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 4
    #>   USUBJID VISITDY PCTPT         NFRLT
    #>   <chr>     <dbl> <chr>         <dbl>
    #> 1 001           1 Pre-dose          0
    #> 2 001           1 1H Post-dose      1
    #> 3 001           1 2H Post-dose      2
    #> 4 001           1 4H Post-dose      4
    #> 5 001           1 24H Post-dose    24

### Single dose study with unit variable

Creating NFRLT with FRLTU unit variable

    derive_var_nfrlt(
      adpc,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT         NFRLT FRLTU
    #>   <chr>     <dbl> <chr>         <dbl> <chr>
    #> 1 001           1 Pre-dose          0 HOURS
    #> 2 001           1 1H Post-dose      1 HOURS
    #> 3 001           1 2H Post-dose      2 HOURS
    #> 4 001           1 4H Post-dose      4 HOURS
    #> 5 001           1 24H Post-dose    24 HOURS

### Single dose study with different output units

Deriving NFRLT in different time units with unit variables

    adpc %>%
      derive_var_nfrlt(
        new_var = NFRLT,
        new_var_unit = FRLTU,
        out_unit = "HOURS",
        tpt_var = PCTPT,
        visit_day = VISITDY
      ) %>%
      derive_var_nfrlt(
        new_var = NFRLTDY,
        new_var_unit = FRLTDYU,
        out_unit = "days",
        tpt_var = PCTPT,
        visit_day = VISITDY
      )
    #> # A tibble: 5 × 7
    #>   USUBJID VISITDY PCTPT         NFRLT FRLTU NFRLTDY FRLTDYU
    #>   <chr>     <dbl> <chr>         <dbl> <chr>   <dbl> <chr>
    #> 1 001           1 Pre-dose          0 HOURS  0      days
    #> 2 001           1 1H Post-dose      1 HOURS  0.0417 days
    #> 3 001           1 2H Post-dose      2 HOURS  0.0833 days
    #> 4 001           1 4H Post-dose      4 HOURS  0.167  days
    #> 5 001           1 24H Post-dose    24 HOURS  1      days   

### Study with screening visits

Handling negative visit days (no Day 0 in clinical trials)

    adpc_screen <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    -14,      "Screening",
      "001",    -7,       "Pre-dose",
      "001",    -1,       "Pre-dose",
      "001",    1,        "Pre-dose",
      "001",    1,        "2H Post-dose"
    )

    derive_var_nfrlt(
      adpc_screen,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT        NFRLT FRLTU
    #>   <chr>     <dbl> <chr>        <dbl> <chr>
    #> 1 001         -14 Screening     -336 HOURS
    #> 2 001          -7 Pre-dose      -168 HOURS
    #> 3 001          -1 Pre-dose       -24 HOURS
    #> 4 001           1 Pre-dose         0 HOURS
    #> 5 001           1 2H Post-dose     2 HOURS

### Multiple dose study

Dosing on Days 1, 8, and 15

    adpc_md <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    1,        "2H Post-dose",
      "001",    8,        "Pre-dose",
      "001",    8,        "2H Post-dose",
      "001",    15,       "Pre-dose",
      "001",    15,       "2H Post-dose"
    )

    derive_var_nfrlt(
      adpc_md,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 6 × 5
    #>   USUBJID VISITDY PCTPT        NFRLT FRLTU
    #>   <chr>     <dbl> <chr>        <dbl> <chr>
    #> 1 001           1 Pre-dose         0 HOURS
    #> 2 001           1 2H Post-dose     2 HOURS
    #> 3 001           8 Pre-dose       168 HOURS
    #> 4 001           8 2H Post-dose   170 HOURS
    #> 5 001          15 Pre-dose       336 HOURS
    #> 6 001          15 2H Post-dose   338 HOURS

### Multiple dose study with days output

Deriving both NFRLT (hours) and NFRLTDY (days) with unit variables

    adpc_md %>%
      derive_var_nfrlt(
        new_var = NFRLT,
        new_var_unit = FRLTU,
        tpt_var = PCTPT,
        visit_day = VISITDY
      ) %>%
      derive_var_nfrlt(
        new_var = NFRLTDY,
        new_var_unit = FRLTDYU,
        out_unit = "days",
        tpt_var = PCTPT,
        visit_day = VISITDY
      )
    #> # A tibble: 6 × 7
    #>   USUBJID VISITDY PCTPT        NFRLT FRLTU NFRLTDY FRLTDYU
    #>   <chr>     <dbl> <chr>        <dbl> <chr>   <dbl> <chr>
    #> 1 001           1 Pre-dose         0 HOURS  0      days
    #> 2 001           1 2H Post-dose     2 HOURS  0.0833 days
    #> 3 001           8 Pre-dose       168 HOURS  7      days
    #> 4 001           8 2H Post-dose   170 HOURS  7.08   days
    #> 5 001          15 Pre-dose       336 HOURS 14      days
    #> 6 001          15 2H Post-dose   338 HOURS 14.1    days   

### Weekly dosing study

Long-term study with weekly dosing, using weeks output

    adpc_weekly <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    8,        "Pre-dose",
      "001",    15,       "Pre-dose",
      "001",    22,       "Pre-dose",
      "001",    29,       "Pre-dose"
    )

    derive_var_nfrlt(
      adpc_weekly,
      new_var = NFRLTWK,
      new_var_unit = FRLTU,
      out_unit = "weeks",
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT    NFRLTWK FRLTU
    #>   <chr>     <dbl> <chr>      <dbl> <chr>
    #> 1 001           1 Pre-dose       0 weeks
    #> 2 001           8 Pre-dose       1 weeks
    #> 3 001          15 Pre-dose       2 weeks
    #> 4 001          22 Pre-dose       3 weeks
    #> 5 001          29 Pre-dose       4 weeks

### Short-term PK study with minutes

Very short timepoints requiring minute precision

    adpc_short <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    1,        "5 MIN POST",
      "001",    1,        "15 MIN POST",
      "001",    1,        "30 MIN POST",
      "001",    1,        "1H POST"
    )

    derive_var_nfrlt(
      adpc_short,
      new_var = NFRLTMIN,
      new_var_unit = FRLTU,
      out_unit = "minutes",
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT       NFRLTMIN FRLTU
    #>   <chr>     <dbl> <chr>          <dbl> <chr>
    #> 1 001           1 Pre-dose           0 minutes
    #> 2 001           1 5 MIN POST         5 minutes
    #> 3 001           1 15 MIN POST       15 minutes
    #> 4 001           1 30 MIN POST       30 minutes
    #> 5 001           1 1H POST           60 minutes

### Custom first dose day

First dose on Day 7 instead of Day 1

    adpc_day7 <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    -1,       "Pre-dose",
      "001",    1,        "Pre-dose",
      "001",    6,        "Pre-dose",
      "001",    7,        "Pre-dose",
      "001",    8,        "Pre-dose"
    )

    derive_var_nfrlt(
      adpc_day7,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      first_dose_day = 7
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT    NFRLT FRLTU
    #>   <chr>     <dbl> <chr>    <dbl> <chr>
    #> 1 001          -1 Pre-dose  -168 HOURS
    #> 2 001           1 Pre-dose  -144 HOURS
    #> 3 001           6 Pre-dose   -24 HOURS
    #> 4 001           7 Pre-dose     0 HOURS
    #> 5 001           8 Pre-dose    24 HOURS

### IV infusion with scalar treatment duration

2-hour infusion duration for all records

    adpc_inf <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    1,        "EOI",
      "001",    1,        "1H Post EOI",
      "001",    1,        "10MIN PRE EOI"
    )

    derive_var_nfrlt(
      adpc_inf,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      treatment_duration = 2
    )
    #> # A tibble: 4 × 5
    #>   USUBJID VISITDY PCTPT         NFRLT FRLTU
    #>   <chr>     <dbl> <chr>         <dbl> <chr>
    #> 1 001           1 Pre-dose       0    HOURS
    #> 2 001           1 EOI            2    HOURS
    #> 3 001           1 1H Post EOI    3    HOURS
    #> 4 001           1 10MIN PRE EOI  1.83 HOURS

### Variable treatment duration

Different treatment durations per subject using a variable

    adpc_var_dur <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,           ~EXDUR,
      "001",    1,        "Pre-dose",       1,
      "001",    1,        "EOI",            1,
      "001",    1,        "1H POST EOI",    1,
      "002",    1,        "Pre-dose",       2,
      "002",    1,        "EOI",            2,
      "002",    1,        "1H POST EOI",    2
    )

    derive_var_nfrlt(
      adpc_var_dur,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      treatment_duration = EXDUR
    )
    #> # A tibble: 6 × 6
    #>   USUBJID VISITDY PCTPT       EXDUR NFRLT FRLTU
    #>   <chr>     <dbl> <chr>       <dbl> <dbl> <chr>
    #> 1 001           1 Pre-dose        1     0 HOURS
    #> 2 001           1 EOI             1     1 HOURS
    #> 3 001           1 1H POST EOI     1     2 HOURS
    #> 4 002           1 Pre-dose        2     0 HOURS
    #> 5 002           1 EOI             2     2 HOURS
    #> 6 002           1 1H POST EOI     2     3 HOURS

### Exposure records without timepoint variable

Deriving NFRLT based only on visit day

    ex <- tribble(
      ~USUBJID, ~VISITDY,
      "001",    1,
      "001",    8,
      "001",    15
    )

    derive_var_nfrlt(
      ex,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      visit_day = VISITDY
    )
    #> # A tibble: 3 × 4
    #>   USUBJID VISITDY NFRLT FRLTU
    #>   <chr>     <dbl> <dbl> <chr>
    #> 1 001           1     0 HOURS
    #> 2 001           8   168 HOURS
    #> 3 001          15   336 HOURS

### Exposure records with different output units

Deriving NFRLT in hours, days, and weeks for exposure records

    ex %>%
      derive_var_nfrlt(
        new_var = NFRLT,
        new_var_unit = FRLTU,
        visit_day = VISITDY
      ) %>%
      derive_var_nfrlt(
        new_var = NFRLTDY,
        new_var_unit = FRLTDYU,
        out_unit = "days",
        visit_day = VISITDY
      ) %>%
      derive_var_nfrlt(
        new_var = NFRLTWK,
        new_var_unit = FRLTWKU,
        out_unit = "weeks",
        visit_day = VISITDY
      )
    #> # A tibble: 3 × 8
    #>   USUBJID VISITDY NFRLT FRLTU NFRLTDY FRLTDYU NFRLTWK FRLTWKU
    #>   <chr>     <dbl> <dbl> <chr>   <dbl> <chr>     <dbl> <chr>
    #> 1 001           1     0 HOURS       0 days          0 weeks
    #> 2 001           8   168 HOURS       7 days          1 weeks
    #> 3 001          15   336 HOURS      14 days          2 weeks  

### Unscheduled visits

Setting NFRLT to NA for unscheduled visits

    adpc_unsched <- tribble(
      ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
      "001",    1,        "VISIT 1",     "Pre-dose",
      "001",    1,        "VISIT 1",     "2H Post-dose",
      "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
      "001",    NA_real_, "UNSCHEDULED", "2H Post-dose"
    )

    derive_var_nfrlt(
      adpc_unsched,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      set_values_to_na = VISIT == "UNSCHEDULED"
    )
    #> # A tibble: 4 × 6
    #>   USUBJID VISITDY VISIT       PCTPT        NFRLT FRLTU
    #>   <chr>     <dbl> <chr>       <chr>        <dbl> <chr>
    #> 1 001           1 VISIT 1     Pre-dose         0 HOURS
    #> 2 001           1 VISIT 1     2H Post-dose     2 HOURS
    #> 3 001          NA UNSCHEDULED Pre-dose        NA <NA>
    #> 4 001          NA UNSCHEDULED 2H Post-dose    NA <NA> 

### Early discontinuation visits

Handling study drug early discontinuation

    adpc_disc <- tribble(
      ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
      "001",    1,        "VISIT 1",                           "Pre-dose",
      "001",    1,        "VISIT 1",                           "2H Post-dose",
      "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "Pre-dose"
    )

    derive_var_nfrlt(
      adpc_disc,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      set_values_to_na = VISIT == "STUDY DRUG EARLY DISCONTINUATION"
    )
    #> # A tibble: 3 × 6
    #>   USUBJID VISITDY VISIT                            PCTPT        NFRLT FRLTU
    #>   <chr>     <dbl> <chr>                            <chr>        <dbl> <chr>
    #> 1 001           1 VISIT 1                          Pre-dose         0 HOURS
    #> 2 001           1 VISIT 1                          2H Post-dose     2 HOURS
    #> 3 001          NA STUDY DRUG EARLY DISCONTINUATION Pre-dose        NA <NA> 

### Multiple exclusion criteria

Excluding multiple visit types

    adpc_multi <- tribble(
      ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
      "001",    1,        "VISIT 1",                           "Pre-dose",
      "001",    NA_real_, "UNSCHEDULED",                       "Pre-dose",
      "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "Pre-dose"
    )

    derive_var_nfrlt(
      adpc_multi,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      set_values_to_na = VISIT %in% c(
        "UNSCHEDULED",
        "STUDY DRUG EARLY DISCONTINUATION"
      )
    )
    #> # A tibble: 3 × 6
    #>   USUBJID VISITDY VISIT                            PCTPT    NFRLT FRLTU
    #>   <chr>     <dbl> <chr>                            <chr>    <dbl> <chr>
    #> 1 001           1 VISIT 1                          Pre-dose     0 HOURS
    #> 2 001          NA UNSCHEDULED                      Pre-dose    NA <NA>
    #> 3 001          NA STUDY DRUG EARLY DISCONTINUATION Pre-dose    NA <NA> 

### Setting special values instead of NA

Using mutate to set NFRLT to 99999 for unscheduled visits

    adpc_unsched_value <- tribble(
      ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
      "001",    1,        "VISIT 1",     "Pre-dose",
      "001",    1,        "VISIT 1",     "2H Post-dose",
      "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
      "001",    NA_real_, "UNSCHEDULED", "2H Post-dose"
    )

    adpc_unsched_value %>%
      derive_var_nfrlt(
        new_var = NFRLT,
        new_var_unit = FRLTU,
        tpt_var = PCTPT,
        visit_day = VISITDY,
        set_values_to_na = VISIT == "UNSCHEDULED"
      ) %>%
      mutate(
        NFRLT = if_else(is.na(NFRLT) & VISIT == "UNSCHEDULED", 99999, NFRLT),
        FRLTU = if_else(is.na(FRLTU) & VISIT == "UNSCHEDULED", "", FRLTU)
      )
    #> # A tibble: 4 × 6
    #>   USUBJID VISITDY VISIT       PCTPT        NFRLT FRLTU
    #>   <chr>     <dbl> <chr>       <chr>        <dbl> <chr>
    #> 1 001           1 VISIT 1     Pre-dose         0 "HOURS"
    #> 2 001           1 VISIT 1     2H Post-dose     2 "HOURS"
    #> 3 001          NA UNSCHEDULED Pre-dose     99999 ""
    #> 4 001          NA UNSCHEDULED 2H Post-dose 99999 ""     

### Custom range method

Using end of range instead of midpoint

    adpc_range <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Pre-dose",
      "001",    1,        "0-6h Post-dose"
    )

    derive_var_nfrlt(
      adpc_range,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      range_method = "end"
    )
    #> # A tibble: 2 × 5
    #>   USUBJID VISITDY PCTPT          NFRLT FRLTU
    #>   <chr>     <dbl> <chr>          <dbl> <chr>
    #> 1 001           1 Pre-dose           0 HOURS
    #> 2 001           1 0-6h Post-dose     6 HOURS

### Alternative terminology

Using "Before" and "After" terminology

    adpc_alt <- tribble(
      ~USUBJID, ~VISITDY, ~PCTPT,
      "001",    1,        "Before",
      "001",    1,        "1H After",
      "001",    1,        "2H After"
    )

    derive_var_nfrlt(
      adpc_alt,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 3 × 5
    #>   USUBJID VISITDY PCTPT    NFRLT FRLTU
    #>   <chr>     <dbl> <chr>    <dbl> <chr>
    #> 1 001           1 Before       0 HOURS
    #> 2 001           1 1H After     1 HOURS
    #> 3 001           1 2H After     2 HOURS

### Reference relative time with RRLTU

Using RRLTU for reference relative time instead of first dose

    derive_var_nfrlt(
      adpc,
      new_var = NRRLT,
      new_var_unit = RRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      first_dose_day = 8
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT         NRRLT RRLTU
    #>   <chr>     <dbl> <chr>         <dbl> <chr>
    #> 1 001           1 Pre-dose       -168 HOURS
    #> 2 001           1 1H Post-dose   -167 HOURS
    #> 3 001           1 2H Post-dose   -166 HOURS
    #> 4 001           1 4H Post-dose   -164 HOURS
    #> 5 001           1 24H Post-dose  -144 HOURS

### Case sensitivity in out_unit

Unit variable preserves the case provided in out_unit

    derive_var_nfrlt(
      adpc,
      new_var = NFRLT,
      new_var_unit = FRLTU,
      out_unit = "HOURS",
      tpt_var = PCTPT,
      visit_day = VISITDY
    )
    #> # A tibble: 5 × 5
    #>   USUBJID VISITDY PCTPT         NFRLT FRLTU
    #>   <chr>     <dbl> <chr>         <dbl> <chr>
    #> 1 001           1 Pre-dose          0 HOURS
    #> 2 001           1 1H Post-dose      1 HOURS
    #> 3 001           1 2H Post-dose      2 HOURS
    #> 4 001           1 4H Post-dose      4 HOURS
    #> 5 001           1 24H Post-dose    24 HOURS
