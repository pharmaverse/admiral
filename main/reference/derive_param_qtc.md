# Adds a Parameter for Corrected QT (an ECG measurement)

Adds a record for corrected QT using either Bazett's, Fridericia's or
Sagie's formula for each by group (e.g., subject and visit) where the
source parameters are available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md).

## Usage

``` r
derive_param_qtc(
  dataset,
  by_vars,
  method,
  set_values_to = default_qtc_paramcd(method),
  qt_code = "QT",
  rr_code = "RR",
  get_unit_expr,
  filter = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `get_unit_expr` arguments
  are expected to be in the dataset. `PARAMCD`, and `AVAL` are expected
  as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset after restricting it by the filter condition
  (`filter` argument) and to the parameters specified by `qt_code` and
  `rr_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  Only variables specified in `by_vars` will be populated in the newly
  created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- method:

  Method used to QT correction

  See
  [`compute_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qtc.md)
  for details.

  Permitted values

  :   `"Bazett"`, `"Fridericia"`, `"Sagie"`

  Default value

  :   none

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations. For example `exprs(PARAMCD = "MAP")` defines the
  parameter code for the new parameter.

  Permitted values

  :   List of variable-value pairs

  Default value

  :   `exprs(PARAMCD = "MAP")`

- qt_code:

  QT parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the QT interval assessments. It is expected that QT is
  measured in ms or msec.

  Permitted values

  :   character value

  Default value

  :   `"QT"`

- rr_code:

  RR parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the RR interval assessments. It is expected that RR is
  measured in ms or msec.

  Permitted values

  :   character value

  Default value

  :   `"RR"`

- get_unit_expr:

  An expression providing the unit of the parameter

  The result is used to check the units of the input parameters.

  Permitted values

  :   An expression which is evaluable in the input dataset and results
      in a character value

  Default value

  :   none

- filter:

  Filter condition

  The specified condition is applied to the input dataset before
  deriving the new parameter, i.e., only observations fulfilling the
  condition are taken into account.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

## Value

The input dataset with the new parameter added. Note, a variable will
only be populated in the new parameter rows if it is specified in
`by_vars`.

## See also

[`compute_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qtc.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

adeg <- tribble(
  ~USUBJID,      ~PARAMCD, ~PARAM,                   ~AVAL, ~AVALU,      ~VISIT,
  "01-701-1015", "HR",     "Heart Rate (beats/min)", 70.14, "beats/min", "BASELINE",
  "01-701-1015", "QT",     "QT Duration (ms)",         370, "ms",        "WEEK 2",
  "01-701-1015", "HR",     "Heart Rate (beats/min)", 62.66, "beats/min", "WEEK 1",
  "01-701-1015", "RR",     "RR Duration (ms)",         710, "ms",        "WEEK 2",
  "01-701-1028", "HR",     "Heart Rate (beats/min)", 85.45, "beats/min", "BASELINE",
  "01-701-1028", "QT",     "QT Duration (ms)",         480, "ms",        "WEEK 2",
  "01-701-1028", "QT",     "QT Duration (ms)",         350, "ms",        "WEEK 3",
  "01-701-1028", "HR",     "Heart Rate (beats/min)", 56.54, "beats/min", "WEEK 3",
  "01-701-1028", "RR",     "RR Duration (ms)",         842, "ms",        "WEEK 2"
)

derive_param_qtc(
  adeg,
  by_vars = exprs(USUBJID, VISIT),
  method = "Bazett",
  set_values_to = exprs(
    PARAMCD = "QTCBR",
    PARAM = "QTcB - Bazett's Correction Formula Rederived (ms)",
    AVALU = "ms"
  ),
  get_unit_expr = AVALU
)
#> # A tibble: 11 × 6
#>    USUBJID     PARAMCD PARAM                                    AVAL AVALU VISIT
#>    <chr>       <chr>   <chr>                                   <dbl> <chr> <chr>
#>  1 01-701-1015 HR      Heart Rate (beats/min)                   70.1 beat… BASE…
#>  2 01-701-1015 QT      QT Duration (ms)                        370   ms    WEEK…
#>  3 01-701-1015 HR      Heart Rate (beats/min)                   62.7 beat… WEEK…
#>  4 01-701-1015 RR      RR Duration (ms)                        710   ms    WEEK…
#>  5 01-701-1028 HR      Heart Rate (beats/min)                   85.4 beat… BASE…
#>  6 01-701-1028 QT      QT Duration (ms)                        480   ms    WEEK…
#>  7 01-701-1028 QT      QT Duration (ms)                        350   ms    WEEK…
#>  8 01-701-1028 HR      Heart Rate (beats/min)                   56.5 beat… WEEK…
#>  9 01-701-1028 RR      RR Duration (ms)                        842   ms    WEEK…
#> 10 01-701-1015 QTCBR   QTcB - Bazett's Correction Formula Red… 439.  ms    WEEK…
#> 11 01-701-1028 QTCBR   QTcB - Bazett's Correction Formula Red… 523.  ms    WEEK…

derive_param_qtc(
  adeg,
  by_vars = exprs(USUBJID, VISIT),
  method = "Fridericia",
  set_values_to = exprs(
    PARAMCD = "QTCFR",
    PARAM = "QTcF - Fridericia's Correction Formula Rederived (ms)",
    AVALU = "ms"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 11 × 6
#>    USUBJID     PARAMCD PARAM                                    AVAL AVALU VISIT
#>    <chr>       <chr>   <chr>                                   <dbl> <chr> <chr>
#>  1 01-701-1015 HR      Heart Rate (beats/min)                   70.1 beat… BASE…
#>  2 01-701-1015 QT      QT Duration (ms)                        370   ms    WEEK…
#>  3 01-701-1015 HR      Heart Rate (beats/min)                   62.7 beat… WEEK…
#>  4 01-701-1015 RR      RR Duration (ms)                        710   ms    WEEK…
#>  5 01-701-1028 HR      Heart Rate (beats/min)                   85.4 beat… BASE…
#>  6 01-701-1028 QT      QT Duration (ms)                        480   ms    WEEK…
#>  7 01-701-1028 QT      QT Duration (ms)                        350   ms    WEEK…
#>  8 01-701-1028 HR      Heart Rate (beats/min)                   56.5 beat… WEEK…
#>  9 01-701-1028 RR      RR Duration (ms)                        842   ms    WEEK…
#> 10 01-701-1015 QTCFR   QTcF - Fridericia's Correction Formula… 415.  ms    WEEK…
#> 11 01-701-1028 QTCFR   QTcF - Fridericia's Correction Formula… 508.  ms    WEEK…

derive_param_qtc(
  adeg,
  by_vars = exprs(USUBJID, VISIT),
  method = "Sagie",
  set_values_to = exprs(
    PARAMCD = "QTLCR",
    PARAM = "QTlc - Sagie's Correction Formula Rederived (ms)",
    AVALU = "ms"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 11 × 6
#>    USUBJID     PARAMCD PARAM                                    AVAL AVALU VISIT
#>    <chr>       <chr>   <chr>                                   <dbl> <chr> <chr>
#>  1 01-701-1015 HR      Heart Rate (beats/min)                   70.1 beat… BASE…
#>  2 01-701-1015 QT      QT Duration (ms)                        370   ms    WEEK…
#>  3 01-701-1015 HR      Heart Rate (beats/min)                   62.7 beat… WEEK…
#>  4 01-701-1015 RR      RR Duration (ms)                        710   ms    WEEK…
#>  5 01-701-1028 HR      Heart Rate (beats/min)                   85.4 beat… BASE…
#>  6 01-701-1028 QT      QT Duration (ms)                        480   ms    WEEK…
#>  7 01-701-1028 QT      QT Duration (ms)                        350   ms    WEEK…
#>  8 01-701-1028 HR      Heart Rate (beats/min)                   56.5 beat… WEEK…
#>  9 01-701-1028 RR      RR Duration (ms)                        842   ms    WEEK…
#> 10 01-701-1015 QTLCR   QTlc - Sagie's Correction Formula Rede… 415.  ms    WEEK…
#> 11 01-701-1028 QTLCR   QTlc - Sagie's Correction Formula Rede… 504.  ms    WEEK…
```
