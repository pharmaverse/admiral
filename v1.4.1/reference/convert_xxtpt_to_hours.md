# Convert `XXTPT` Strings to Hours

**\[experimental\]**

Converts CDISC timepoint strings (e.g., `PCTPT`, `VSTPT`, `EGTPT`,
`ISTPT`, `LBTPT`) into numeric hours for analysis. The function handles
common dose-centric formats including pre-dose, post-dose
(hours/minutes), days, time ranges, and treatment-related time markers.

## Usage

``` r
convert_xxtpt_to_hours(
  xxtpt,
  treatment_duration = 0,
  range_method = "midpoint"
)
```

## Arguments

- xxtpt:

  A character vector of timepoint descriptions from SDTM `--TPT`
  variables (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`). Can
  contain `NA` values.

  Permitted values

  :   character vector

  Default value

  :   none

- treatment_duration:

  Numeric value(s) specifying the duration of treatment in hours. Used
  to convert "EOI/EOT" (End of Infusion/Treatment) patterns and patterns
  describing time after end of treatment. Must be non-negative. Can be
  either:

  - A single value (used for all timepoints), or

  - A vector of the same length as `xxtpt` (one value per timepoint)

  Default is 0 hours (for instantaneous treatments like oral
  medications).

  Permitted values

  :   numeric scalar or vector (non-negative)

  Default value

  :   `0`

- range_method:

  Method for converting time ranges to single values. Options are
  "midpoint" (default), "start", or "end". For example, "0-6h" with
  midpoint returns 3, with start returns 0, with end returns 6.

  Permitted values

  :   character scalar ("midpoint", "start", or "end")

  Default value

  :   `"midpoint"`

## Value

A numeric vector of timepoints in hours. Returns `NA_real_` for:

- Input `NA` values

- Unrecognized timepoint formats

- Non-time descriptors (e.g., "Morning", "Evening")

Returns `numeric(0)` for empty input.

## Details

The function recognizes the following patterns (all case-insensitive):

**Special Cases:**

- `"Screening"` -\> 0

- `"Pre-dose"`, `"Predose"`, `"Pre-treatment"`, `"Pre-infusion"`,
  `"Pre-inf"`, `"Before"`, `"Infusion"`, `"0H"` -\> 0

- `"EOI"`, `"EOT"`, `"End of Infusion"`, `"End of Treatment"`,
  `"After End of Infusion"`, `"After End of Treatment"` -\>
  `treatment_duration` (default: 0)

- `"Morning"`, `"Evening"` -\> `NA_real_`

- Unrecognized values -\> `NA_real_`

**Time Ranges:** Time ranges are converted based on the `range_method`
parameter:

- `"0-6h Post-dose"` with `range_method = "midpoint"` (default) -\> 3

- `"0-6h Post-dose"` with `range_method = "start"` -\> 0

- `"0-6h Post-dose"` with `range_method = "end"` -\> 6

- `"0-4H PRIOR START OF INFUSION"` with midpoint -\> -2 (negative for
  prior)

- `"8-16H POST START OF INFUSION"` with midpoint -\> 12

- `"0-4H AFTER EOI"` with midpoint and treatment_duration=1 -\> 3 (1 +
  2)

- `"0-4H EOT"` with midpoint and treatment_duration=0 -\> 2

- `"4-8H AFTER END OF INFUSION"` with midpoint and treatment_duration=1
  -\> 7 (1 + 6)

- `"4-8H POST INFUSION"` with midpoint and treatment_duration=1 -\> 7
  (1 + 6)

- `"4-8H POST-INF"` with midpoint and treatment_duration=1 -\> 7 (1 + 6)

**Time-based Conversions:**

- **Days**: `"Day 1"` -\> 24, `"2D"` -\> 48, `"30 DAYS AFTER LAST"` -\>
  720 (requires unit indicator; bare numbers like `"2"` return `NA`)

- **Hours + Minutes**: `"1H30M"` -\> 1.5

- **Hours**: `"2 hours"` -\> 2, `"1 HOUR POST"` -\> 1

- **Minutes**: `"30M"` -\> 0.5, `"30 MIN POST"` -\> 0.5

- **Predose**: `"5 MIN PREDOSE"` -\> -0.0833, `"5 MIN PRE-DOSE"` -\>
  -0.0833

- **Before treatment**: `"5 MIN BEFORE"` -\> -0.0833

- **Post EOI/EOT**: `"1 HOUR POST EOI"` -\> treatment_duration + 1,
  `"24 HR POST INF"` -\> treatment_duration + 24, `"24 HR POST-INF"` -\>
  treatment_duration + 24, `"1 HOUR AFTER EOT"` -\> treatment_duration +
  1

- **After end**: `"30MIN AFTER END OF INFUSION"` -\>
  treatment_duration + 0.5

- **Start of infusion/treatment**: `"8H PRIOR START OF INFUSION"` -\>
  -8, `"8H BEFORE START OF TREATMENT"` -\> -8

- **Pre EOI/EOT**: `"10MIN PRE EOI"` -\> treatment_duration - 1/6,
  `"10MIN BEFORE EOT"` -\> treatment_duration - 1/6

**Supported Unit Formats:**

- Hours: H, h, HR, hr, HOUR, hour (with optional plurals)

- Minutes: M, m, MIN, min, MINUTE, minute (with optional plurals)

- Days: D, d, DAY, day (with optional plurals)

- Flexible whitespace and optional "Post-dose", "POST", "After last"
  suffixes

- Hyphens in compound terms: "PRE-DOSE", "POST-INF", "POST-INFUSION"

**Understanding POST/AFTER Patterns:**

It's important to distinguish between patterns relative to treatment
**start** versus treatment **end**:

- **Relative to START** (treatment_duration NOT added):

  - `"1H POST"`, `"1H AFTER"`, `"30M POST"` -\> Time from dose/treatment
    start

  - These patterns assume treatment starts at time 0

  - Example: `"1H POST"` -\> 1 hour (regardless of treatment_duration)

- **Relative to END** (treatment_duration IS added):

  - `"1H POST EOI"`, `"1H AFTER EOT"`, `"1H POST INFUSION"` -\> Time
    from treatment end

  - These patterns account for when treatment ends (start + duration)

  - Example: `"1H POST EOI"` with treatment_duration=2 -\> 3 hours (2 +
    1)

This distinction follows standard pharmacokinetic conventions where
"post-dose" refers to time from treatment initiation, while "post end of
infusion" refers to time from treatment completion.

**Vectorized Treatment Duration:**

When `treatment_duration` is a vector, each timepoint uses its
corresponding treatment duration value. This is useful when different
records have different treatment durations (e.g., different infusion
lengths).

## See also

Date/Time Computation Functions that returns a vector:
[`compute_age_years()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_age_years.md),
[`compute_dtf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_dtf.md),
[`compute_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_duration.md),
[`compute_tmf()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_tmf.md),
[`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_date_to_dtm.md),
[`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dt.md),
[`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_dtc_to_dtm.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dt.md),
[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dtm.md)

## Examples

### Basic timepoint patterns

Convert basic dose-centric patterns to hours

    convert_xxtpt_to_hours(c(
      "Screening",
      "Pre-dose",
      "Pre-treatment",
      "Before",
      "30M",
      "1H",
      "2H POSTDOSE",
      "Day 1"
    ))
    #> [1]  0.0  0.0  0.0  0.0  0.5  1.0  2.0 24.0

### Predose and before patterns

Convert predose/before patterns that return negative times

    convert_xxtpt_to_hours(c("5 MIN PREDOSE", "5 MIN PRE-DOSE", "1 HOUR BEFORE"))
    #> [1] -0.08333333 -0.08333333 -1.00000000

### Treatment-related patterns (oral medications)

With default treatment_duration = 0 for oral medications

    convert_xxtpt_to_hours(c(
      "EOT",
      "1 HOUR POST EOT",
      "1 HOUR AFTER EOT",
      "After End of Treatment"
    ))
    #> [1] 0 1 1 0

### Infusion-related patterns

With treatment_duration = 1 hour for IV infusions

    convert_xxtpt_to_hours(
      c(
        "EOI",
        "1 HOUR POST EOI",
        "24 HR POST INF",
        "24 HR POST-INF",
        "30MIN AFTER END OF INFUSION",
        "8H PRIOR START OF INFUSION",
        "10MIN PRE EOI"
      ),
      treatment_duration = 1
    )
    #> [1]  1.0000000  2.0000000 25.0000000 25.0000000  1.5000000 -8.0000000  0.8333333

### Vectorized treatment duration

Different treatment durations per timepoint

    convert_xxtpt_to_hours(
      c("EOI", "1 HOUR POST EOI", "EOI", "1 HOUR POST EOI"),
      treatment_duration = c(1, 1, 2, 2)
    )
    #> [1] 1 2 2 3

### Time ranges with midpoint method

Default midpoint method for ranges

    convert_xxtpt_to_hours(c(
      "0-6h Post-dose",
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION"
    ))
    #> [1]  3 -2 12

### Time ranges with custom methods

Specify start or end method for ranges

    convert_xxtpt_to_hours("0-6h Post-dose", range_method = "end")
    #> [1] 6
    convert_xxtpt_to_hours("0-6h Post-dose", range_method = "start")
    #> [1] 0

### Ranges relative to EOI/EOT

Time ranges after end of infusion/treatment

    convert_xxtpt_to_hours(
      c(
        "0-4H AFTER EOI",
        "0-4H POST EOI",
        "4-8H AFTER END OF INFUSION",
        "4-8H AFTER EOT",
        "4-8H POST INFUSION",
        "4-8H POST-INF"
      ),
      treatment_duration = 1
    )
    #> [1] 3 3 7 7 7 7

### POST vs POST EOI distinction

Difference between POST (from start) and POST EOI (from end)

    convert_xxtpt_to_hours(
      c("Pre-dose", "1H POST", "2H POST", "4H POST"),
      treatment_duration = 2
    )
    #> [1] 0 1 2 4

    convert_xxtpt_to_hours(
      c("Pre-dose", "EOI", "1H POST EOI", "2H POST EOI"),
      treatment_duration = 2
    )
    #> [1] 0 2 3 4

    convert_xxtpt_to_hours(
      c("1H POST", "1H POST EOI", "1H POST INFUSION"),
      treatment_duration = 2
    )
    #> [1] 1 3 3
