# Adds a Parameter for Mean Arterial Pressure

Adds a record for mean arterial pressure (MAP) for each by group (e.g.,
subject and visit) where the source parameters are available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_computed.md).

## Usage

``` r
derive_param_map(
  dataset,
  by_vars,
  set_values_to = exprs(PARAMCD = "MAP"),
  sysbp_code = "SYSBP",
  diabp_code = "DIABP",
  hr_code = NULL,
  get_unit_expr,
  filter = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD`, and `AVAL` are expected as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset after restricting it by the filter condition
  (`filter` parameter) and to the parameters specified by `sysbp_code`,
  `diabp_code` and `hr_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` an observation is added to the
  output dataset. Only variables specified in `by_vars` will be
  populated in the newly created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

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

- sysbp_code:

  Systolic blood pressure parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the systolic blood pressure assessments.

  Permitted values

  :   character value

  Default value

  :   `"SYSBP"`

- diabp_code:

  Diastolic blood pressure parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the diastolic blood pressure assessments.

  Permitted values

  :   character value

  Default value

  :   `"DIABP"`

- hr_code:

  Heart rate parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the heart rate assessments.

  Permitted values

  :   character value

  Default value

  :   `NULL`

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

## Details

The analysis value of the new parameter is derived as \$\$\frac{2DIABP +
SYSBP}{3}\$\$ if it is based on diastolic and systolic blood pressure
and \$\$DIABP + 0.01 e^{4.14 - \frac{40.74}{HR}} (SYSBP - DIABP)\$\$ if
it is based on diastolic, systolic blood pressure, and heart rate.

## See also

[`compute_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_map.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_framingham.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
  "01-701-1015", "PULSE", "Pulse (beats/min)", 59, "BASELINE",
  "01-701-1015", "PULSE", "Pulse (beats/min)", 61, "WEEK 2",
  "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
  "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "WEEK 2",
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "WEEK 2",
  "01-701-1028", "PULSE", "Pulse (beats/min)", 62, "BASELINE",
  "01-701-1028", "PULSE", "Pulse (beats/min)", 77, "WEEK 2",
  "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "BASELINE",
  "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "WEEK 2",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "BASELINE",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "WEEK 2"
)

# Derive MAP based on diastolic and systolic blood pressure
advs %>%
  derive_param_map(
    by_vars = exprs(USUBJID, VISIT),
    set_values_to = exprs(
      PARAMCD = "MAP",
      PARAM = "Mean Arterial Pressure (mmHg)"
    ),
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  filter(PARAMCD != "PULSE")
#> # A tibble: 12 × 5
#>    USUBJID     PARAMCD PARAM                            AVAL VISIT   
#>    <chr>       <chr>   <chr>                           <dbl> <chr>   
#>  1 01-701-1015 DIABP   Diastolic Blood Pressure (mmHg)  51   BASELINE
#>  2 01-701-1015 DIABP   Diastolic Blood Pressure (mmHg)  50   WEEK 2  
#>  3 01-701-1015 SYSBP   Systolic Blood Pressure (mmHg)  121   BASELINE
#>  4 01-701-1015 SYSBP   Systolic Blood Pressure (mmHg)  121   WEEK 2  
#>  5 01-701-1028 DIABP   Diastolic Blood Pressure (mmHg)  79   BASELINE
#>  6 01-701-1028 DIABP   Diastolic Blood Pressure (mmHg)  80   WEEK 2  
#>  7 01-701-1028 SYSBP   Systolic Blood Pressure (mmHg)  130   BASELINE
#>  8 01-701-1028 SYSBP   Systolic Blood Pressure (mmHg)  132   WEEK 2  
#>  9 01-701-1015 MAP     Mean Arterial Pressure (mmHg)    74.3 BASELINE
#> 10 01-701-1015 MAP     Mean Arterial Pressure (mmHg)    73.7 WEEK 2  
#> 11 01-701-1028 MAP     Mean Arterial Pressure (mmHg)    96   BASELINE
#> 12 01-701-1028 MAP     Mean Arterial Pressure (mmHg)    97.3 WEEK 2  

# Derive MAP based on diastolic and systolic blood pressure and heart rate
derive_param_map(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  hr_code = "PULSE",
  set_values_to = exprs(
    PARAMCD = "MAP",
    PARAM = "Mean Arterial Pressure (mmHg)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 16 × 5
#>    USUBJID     PARAMCD PARAM                            AVAL VISIT   
#>    <chr>       <chr>   <chr>                           <dbl> <chr>   
#>  1 01-701-1015 PULSE   Pulse (beats/min)                59   BASELINE
#>  2 01-701-1015 PULSE   Pulse (beats/min)                61   WEEK 2  
#>  3 01-701-1015 DIABP   Diastolic Blood Pressure (mmHg)  51   BASELINE
#>  4 01-701-1015 DIABP   Diastolic Blood Pressure (mmHg)  50   WEEK 2  
#>  5 01-701-1015 SYSBP   Systolic Blood Pressure (mmHg)  121   BASELINE
#>  6 01-701-1015 SYSBP   Systolic Blood Pressure (mmHg)  121   WEEK 2  
#>  7 01-701-1028 PULSE   Pulse (beats/min)                62   BASELINE
#>  8 01-701-1028 PULSE   Pulse (beats/min)                77   WEEK 2  
#>  9 01-701-1028 DIABP   Diastolic Blood Pressure (mmHg)  79   BASELINE
#> 10 01-701-1028 DIABP   Diastolic Blood Pressure (mmHg)  80   WEEK 2  
#> 11 01-701-1028 SYSBP   Systolic Blood Pressure (mmHg)  130   BASELINE
#> 12 01-701-1028 SYSBP   Systolic Blood Pressure (mmHg)  132   WEEK 2  
#> 13 01-701-1015 MAP     Mean Arterial Pressure (mmHg)    73.0 BASELINE
#> 14 01-701-1015 MAP     Mean Arterial Pressure (mmHg)    72.9 WEEK 2  
#> 15 01-701-1028 MAP     Mean Arterial Pressure (mmHg)    95.6 BASELINE
#> 16 01-701-1028 MAP     Mean Arterial Pressure (mmHg)    99.2 WEEK 2  
```
