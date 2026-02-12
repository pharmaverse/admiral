# Adds a Parameter for Derived RR (an ECG measurement)

Adds a record for derived RR based on heart rate for each by group
(e.g., subject and visit) where the source parameters are available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md).

The analysis value of the new parameter is derived as
\$\$\frac{60000}{HR}\$\$

## Usage

``` r
derive_param_rr(
  dataset,
  by_vars,
  set_values_to = exprs(PARAMCD = "RRR"),
  hr_code = "HR",
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
  (`filter` argument) and to the parameters specified by `hr_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` an observation is added to the
  output dataset. Only variables specified in `by_vars` will be
  populated in the newly created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
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

- hr_code:

  HR parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the heart rate assessments.

  Permitted values

  :   character value

  Default value

  :   `"HR"`

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

[`compute_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_rr.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_qtc.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

adeg <- tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
  "01-701-1015", "HR", "Heart Rate", 70.14, "beats/min", "BASELINE",
  "01-701-1015", "QT", "QT Duration", 370, "ms", "WEEK 2",
  "01-701-1015", "HR", "Heart Rate", 62.66, "beats/min", "WEEK 1",
  "01-701-1015", "RR", "RR Duration", 710, "ms", "WEEK 2",
  "01-701-1028", "HR", "Heart Rate", 85.45, "beats/min", "BASELINE",
  "01-701-1028", "QT", "QT Duration", 480, "ms", "WEEK 2",
  "01-701-1028", "QT", "QT Duration", 350, "ms", "WEEK 3",
  "01-701-1028", "HR", "Heart Rate", 56.54, "beats/min", "WEEK 3",
  "01-701-1028", "RR", "RR Duration", 842, "ms", "WEEK 2"
)

derive_param_rr(
  adeg,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(
    PARAMCD = "RRR",
    PARAM = "RR Duration Rederived (ms)",
    AVALU = "ms"
  ),
  get_unit_expr = AVALU
)
#> # A tibble: 13 × 6
#>    USUBJID     PARAMCD PARAM                        AVAL AVALU     VISIT   
#>    <chr>       <chr>   <chr>                       <dbl> <chr>     <chr>   
#>  1 01-701-1015 HR      Heart Rate                   70.1 beats/min BASELINE
#>  2 01-701-1015 QT      QT Duration                 370   ms        WEEK 2  
#>  3 01-701-1015 HR      Heart Rate                   62.7 beats/min WEEK 1  
#>  4 01-701-1015 RR      RR Duration                 710   ms        WEEK 2  
#>  5 01-701-1028 HR      Heart Rate                   85.4 beats/min BASELINE
#>  6 01-701-1028 QT      QT Duration                 480   ms        WEEK 2  
#>  7 01-701-1028 QT      QT Duration                 350   ms        WEEK 3  
#>  8 01-701-1028 HR      Heart Rate                   56.5 beats/min WEEK 3  
#>  9 01-701-1028 RR      RR Duration                 842   ms        WEEK 2  
#> 10 01-701-1015 RRR     RR Duration Rederived (ms)  855.  ms        BASELINE
#> 11 01-701-1015 RRR     RR Duration Rederived (ms)  958.  ms        WEEK 1  
#> 12 01-701-1028 RRR     RR Duration Rederived (ms)  702.  ms        BASELINE
#> 13 01-701-1028 RRR     RR Duration Rederived (ms) 1061.  ms        WEEK 3  
```
