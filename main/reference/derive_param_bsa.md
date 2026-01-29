# Adds a Parameter for BSA (Body Surface Area) Using the Specified Method

Adds a record for BSA (Body Surface Area) using the specified derivation
method for each by group (e.g., subject and visit) where the source
parameters are available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md).

## Usage

``` r
derive_param_bsa(
  dataset,
  by_vars,
  method,
  set_values_to = exprs(PARAMCD = "BSA"),
  height_code = "HEIGHT",
  weight_code = "WEIGHT",
  get_unit_expr,
  filter = NULL,
  constant_by_vars = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD`, and `AVAL` are expected as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset after restricting it by the filter condition
  (`filter` parameter) and to the parameters specified by `HEIGHT` and
  `WEIGHT`.

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` an observation is added to the
  output dataset. Only variables specified in `by_vars` will be
  populated in the newly created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- method:

  Derivation method to use. Note that `HEIGHT` is expected in cm and
  `WEIGHT` is expected in kg:

  Mosteller: `sqrt(height * weight / 3600)`

  DuBois-DuBois: `0.20247 * (height/100) ^ 0.725 * weight ^ 0.425`

  Haycock: `0.024265 * height ^ 0.3964 * weight ^ 0.5378`

  Gehan-George: `0.0235 * height ^ 0.42246 * weight ^ 0.51456`

  Boyd:
  `0.0003207 * (height ^ 0.3) * (1000 * weight) ^ (0.7285 - (0.0188 * log10(1000 * weight)))`

  Fujimoto: `0.008883 * height ^ 0.663 * weight ^ 0.444`

  Takahira: `0.007241 * height ^ 0.725 * weight ^ 0.425`

  Permitted values

  :   character value

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

- height_code:

  HEIGHT parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the HEIGHT assessments. It is expected that HEIGHT is
  measured in cm.

  Permitted values

  :   character value

  Default value

  :   `"HEIGHT"`

- weight_code:

  WEIGHT parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the WEIGHT assessments. It is expected that WEIGHT is
  measured in kg.

  Permitted values

  :   character value

  Default value

  :   `"WEIGHT"`

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

- constant_by_vars:

  By variables for when HEIGHT is constant

  When HEIGHT is constant, the HEIGHT parameters (measured only once)
  are merged to the other parameters using the specified variables.

  If height is constant (e.g. only measured once at screening or
  baseline) then use `constant_by_vars` to select the subject-level
  variable to merge on (e.g. `USUBJID`). This will produce BSA at all
  visits where weight is measured. Otherwise it will only be calculated
  at visits with both height and weight collected.

  Default value

  :   `NULL`

## Value

The input dataset with the new parameter added. Note, a variable will
only be populated in the new parameter rows if it is specified in
`by_vars`.

## See also

[`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

# Example 1: Derive BSA where height is measured only once using constant_by_vars
advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
  "01-701-1015", "HEIGHT", "Height (cm)", 170, "BASELINE",
  "01-701-1015", "WEIGHT", "Weight (kg)", 75, "BASELINE",
  "01-701-1015", "WEIGHT", "Weight (kg)", 78, "MONTH 1",
  "01-701-1015", "WEIGHT", "Weight (kg)", 80, "MONTH 2",
  "01-701-1028", "HEIGHT", "Height (cm)", 185, "BASELINE",
  "01-701-1028", "WEIGHT", "Weight (kg)", 90, "BASELINE",
  "01-701-1028", "WEIGHT", "Weight (kg)", 88, "MONTH 1",
  "01-701-1028", "WEIGHT", "Weight (kg)", 85, "MONTH 2",
)

derive_param_bsa(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  method = "Mosteller",
  set_values_to = exprs(
    PARAMCD = "BSA",
    PARAM = "Body Surface Area (m^2)"
  ),
  get_unit_expr = extract_unit(PARAM),
  constant_by_vars = exprs(USUBJID)
)
#> # A tibble: 14 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL VISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>   
#>  1 01-701-1015 HEIGHT  Height (cm)             170    BASELINE
#>  2 01-701-1015 WEIGHT  Weight (kg)              75    BASELINE
#>  3 01-701-1015 WEIGHT  Weight (kg)              78    MONTH 1 
#>  4 01-701-1015 WEIGHT  Weight (kg)              80    MONTH 2 
#>  5 01-701-1028 HEIGHT  Height (cm)             185    BASELINE
#>  6 01-701-1028 WEIGHT  Weight (kg)              90    BASELINE
#>  7 01-701-1028 WEIGHT  Weight (kg)              88    MONTH 1 
#>  8 01-701-1028 WEIGHT  Weight (kg)              85    MONTH 2 
#>  9 01-701-1015 BSA     Body Surface Area (m^2)   1.88 BASELINE
#> 10 01-701-1015 BSA     Body Surface Area (m^2)   1.92 MONTH 1 
#> 11 01-701-1015 BSA     Body Surface Area (m^2)   1.94 MONTH 2 
#> 12 01-701-1028 BSA     Body Surface Area (m^2)   2.15 BASELINE
#> 13 01-701-1028 BSA     Body Surface Area (m^2)   2.13 MONTH 1 
#> 14 01-701-1028 BSA     Body Surface Area (m^2)   2.09 MONTH 2 

derive_param_bsa(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  method = "Fujimoto",
  set_values_to = exprs(
    PARAMCD = "BSA",
    PARAM = "Body Surface Area (m^2)"
  ),
  get_unit_expr = extract_unit(PARAM),
  constant_by_vars = exprs(USUBJID)
)
#> # A tibble: 14 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL VISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>   
#>  1 01-701-1015 HEIGHT  Height (cm)             170    BASELINE
#>  2 01-701-1015 WEIGHT  Weight (kg)              75    BASELINE
#>  3 01-701-1015 WEIGHT  Weight (kg)              78    MONTH 1 
#>  4 01-701-1015 WEIGHT  Weight (kg)              80    MONTH 2 
#>  5 01-701-1028 HEIGHT  Height (cm)             185    BASELINE
#>  6 01-701-1028 WEIGHT  Weight (kg)              90    BASELINE
#>  7 01-701-1028 WEIGHT  Weight (kg)              88    MONTH 1 
#>  8 01-701-1028 WEIGHT  Weight (kg)              85    MONTH 2 
#>  9 01-701-1015 BSA     Body Surface Area (m^2)   1.82 BASELINE
#> 10 01-701-1015 BSA     Body Surface Area (m^2)   1.85 MONTH 1 
#> 11 01-701-1015 BSA     Body Surface Area (m^2)   1.87 MONTH 2 
#> 12 01-701-1028 BSA     Body Surface Area (m^2)   2.09 BASELINE
#> 13 01-701-1028 BSA     Body Surface Area (m^2)   2.07 MONTH 1 
#> 14 01-701-1028 BSA     Body Surface Area (m^2)   2.03 MONTH 2 

# Example 2: Derive BSA where height is measured only once and keep only one record
# where both height and weight are measured.

derive_param_bsa(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  method = "Mosteller",
  set_values_to = exprs(
    PARAMCD = "BSA",
    PARAM = "Body Surface Area (m^2)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 10 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL VISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>   
#>  1 01-701-1015 HEIGHT  Height (cm)             170    BASELINE
#>  2 01-701-1015 WEIGHT  Weight (kg)              75    BASELINE
#>  3 01-701-1015 WEIGHT  Weight (kg)              78    MONTH 1 
#>  4 01-701-1015 WEIGHT  Weight (kg)              80    MONTH 2 
#>  5 01-701-1028 HEIGHT  Height (cm)             185    BASELINE
#>  6 01-701-1028 WEIGHT  Weight (kg)              90    BASELINE
#>  7 01-701-1028 WEIGHT  Weight (kg)              88    MONTH 1 
#>  8 01-701-1028 WEIGHT  Weight (kg)              85    MONTH 2 
#>  9 01-701-1015 BSA     Body Surface Area (m^2)   1.88 BASELINE
#> 10 01-701-1028 BSA     Body Surface Area (m^2)   2.15 BASELINE

# Example 3: Pediatric study where height and weight are measured multiple times
advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
  "01-101-1001", "HEIGHT", "Height (cm)", 47.1, "BASELINE",
  "01-101-1001", "HEIGHT", "Height (cm)", 59.1, "WEEK 12",
  "01-101-1001", "HEIGHT", "Height (cm)", 64.7, "WEEK 24",
  "01-101-1001", "HEIGHT", "Height (cm)", 68.2, "WEEK 48",
  "01-101-1001", "WEIGHT", "Weight (kg)", 2.6, "BASELINE",
  "01-101-1001", "WEIGHT", "Weight (kg)", 5.3, "WEEK 12",
  "01-101-1001", "WEIGHT", "Weight (kg)", 6.7, "WEEK 24",
  "01-101-1001", "WEIGHT", "Weight (kg)", 7.4, "WEEK 48",
)
derive_param_bsa(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  method = "Mosteller",
  set_values_to = exprs(
    PARAMCD = "BSA",
    PARAM = "Body Surface Area (m^2)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 12 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL VISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>   
#>  1 01-101-1001 HEIGHT  Height (cm)             47.1   BASELINE
#>  2 01-101-1001 HEIGHT  Height (cm)             59.1   WEEK 12 
#>  3 01-101-1001 HEIGHT  Height (cm)             64.7   WEEK 24 
#>  4 01-101-1001 HEIGHT  Height (cm)             68.2   WEEK 48 
#>  5 01-101-1001 WEIGHT  Weight (kg)              2.6   BASELINE
#>  6 01-101-1001 WEIGHT  Weight (kg)              5.3   WEEK 12 
#>  7 01-101-1001 WEIGHT  Weight (kg)              6.7   WEEK 24 
#>  8 01-101-1001 WEIGHT  Weight (kg)              7.4   WEEK 48 
#>  9 01-101-1001 BSA     Body Surface Area (m^2)  0.184 BASELINE
#> 10 01-101-1001 BSA     Body Surface Area (m^2)  0.295 WEEK 12 
#> 11 01-101-1001 BSA     Body Surface Area (m^2)  0.347 WEEK 24 
#> 12 01-101-1001 BSA     Body Surface Area (m^2)  0.374 WEEK 48 
```
