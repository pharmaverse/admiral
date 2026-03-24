# Adds a Parameter for BMI

Adds a record for BMI/Body Mass Index using Weight and Height each by
group (e.g., subject and visit) where the source parameters are
available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md).

## Usage

``` r
derive_param_bmi(
  dataset,
  by_vars,
  set_values_to = exprs(PARAMCD = "BMI"),
  weight_code = "WEIGHT",
  height_code = "HEIGHT",
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
  (`filter` parameter) and to the parameters specified by `weight_code`
  and `height_code`.

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

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations. For example `exprs(PARAMCD = "MAP")` defines the
  parameter code for the new parameter.

  Permitted values

  :   List of variable-value pairs

  Default value

  :   `exprs(PARAMCD = "MAP")`

- weight_code:

  WEIGHT parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the WEIGHT. It is expected that WEIGHT is measured in kg

  Permitted values

  :   character value

  Default value

  :   `"WEIGHT"`

- height_code:

  HEIGHT parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the HEIGHT. It is expected that HEIGHT is measured in cm

  Permitted values

  :   logical scalar

  Default value

  :   `"HEIGHT"`

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
  variable to merge on (e.g. `USUBJID`). This will produce BMI at all
  visits where weight is measured. Otherwise it will only be calculated
  at visits with both height and weight collected.

  Default value

  :   `NULL`

## Value

The input dataset with the new parameter added. Note, a variable will
only be populated in the new parameter rows if it is specified in
`by_vars`.

## Details

The analysis value of the new parameter is derived as \$\$BMI =
\frac{WEIGHT}{HEIGHT^2}\$\$

## See also

[`compute_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bmi.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
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
# Example 1: Derive BMI where height is measured only once using constant_by_vars
advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISIT,
  "01-701-1015", "HEIGHT", "Height (cm)", 147, "SCREENING",
  "01-701-1015", "WEIGHT", "Weight (kg)", 54.0, "SCREENING",
  "01-701-1015", "WEIGHT", "Weight (kg)", 54.4, "BASELINE",
  "01-701-1015", "WEIGHT", "Weight (kg)", 53.1, "WEEK 2",
  "01-701-1028", "HEIGHT", "Height (cm)", 163, "SCREENING",
  "01-701-1028", "WEIGHT", "Weight (kg)", 78.5, "SCREENING",
  "01-701-1028", "WEIGHT", "Weight (kg)", 80.3, "BASELINE",
  "01-701-1028", "WEIGHT", "Weight (kg)", 80.7, "WEEK 2"
)

derive_param_bmi(
  advs,
  by_vars = exprs(USUBJID, AVISIT),
  weight_code = "WEIGHT",
  height_code = "HEIGHT",
  set_values_to = exprs(
    PARAMCD = "BMI",
    PARAM = "Body Mass Index (kg/m^2)"
  ),
  get_unit_expr = extract_unit(PARAM),
  constant_by_vars = exprs(USUBJID)
)
#> # A tibble: 14 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL AVISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>    
#>  1 01-701-1015 HEIGHT  Height (cm)              147   SCREENING
#>  2 01-701-1015 WEIGHT  Weight (kg)               54   SCREENING
#>  3 01-701-1015 WEIGHT  Weight (kg)               54.4 BASELINE 
#>  4 01-701-1015 WEIGHT  Weight (kg)               53.1 WEEK 2   
#>  5 01-701-1028 HEIGHT  Height (cm)              163   SCREENING
#>  6 01-701-1028 WEIGHT  Weight (kg)               78.5 SCREENING
#>  7 01-701-1028 WEIGHT  Weight (kg)               80.3 BASELINE 
#>  8 01-701-1028 WEIGHT  Weight (kg)               80.7 WEEK 2   
#>  9 01-701-1015 BMI     Body Mass Index (kg/m^2)  25.0 SCREENING
#> 10 01-701-1015 BMI     Body Mass Index (kg/m^2)  25.2 BASELINE 
#> 11 01-701-1015 BMI     Body Mass Index (kg/m^2)  24.6 WEEK 2   
#> 12 01-701-1028 BMI     Body Mass Index (kg/m^2)  29.5 SCREENING
#> 13 01-701-1028 BMI     Body Mass Index (kg/m^2)  30.2 BASELINE 
#> 14 01-701-1028 BMI     Body Mass Index (kg/m^2)  30.4 WEEK 2   

# Example 2: Derive BMI where height is measured only once and keep only one record
# where both height and weight are measured.
derive_param_bmi(
  advs,
  by_vars = exprs(USUBJID, AVISIT),
  weight_code = "WEIGHT",
  height_code = "HEIGHT",
  set_values_to = exprs(
    PARAMCD = "BMI",
    PARAM = "Body Mass Index (kg/m^2)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 10 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL AVISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>    
#>  1 01-701-1015 HEIGHT  Height (cm)              147   SCREENING
#>  2 01-701-1015 WEIGHT  Weight (kg)               54   SCREENING
#>  3 01-701-1015 WEIGHT  Weight (kg)               54.4 BASELINE 
#>  4 01-701-1015 WEIGHT  Weight (kg)               53.1 WEEK 2   
#>  5 01-701-1028 HEIGHT  Height (cm)              163   SCREENING
#>  6 01-701-1028 WEIGHT  Weight (kg)               78.5 SCREENING
#>  7 01-701-1028 WEIGHT  Weight (kg)               80.3 BASELINE 
#>  8 01-701-1028 WEIGHT  Weight (kg)               80.7 WEEK 2   
#>  9 01-701-1015 BMI     Body Mass Index (kg/m^2)  25.0 SCREENING
#> 10 01-701-1028 BMI     Body Mass Index (kg/m^2)  29.5 SCREENING

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

derive_param_bmi(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  weight_code = "WEIGHT",
  height_code = "HEIGHT",
  set_values_to = exprs(
    PARAMCD = "BMI",
    PARAM = "Body Mass Index (kg/m^2)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 12 × 5
#>    USUBJID     PARAMCD PARAM                     AVAL VISIT   
#>    <chr>       <chr>   <chr>                    <dbl> <chr>   
#>  1 01-101-1001 HEIGHT  Height (cm)               47.1 BASELINE
#>  2 01-101-1001 HEIGHT  Height (cm)               59.1 WEEK 12 
#>  3 01-101-1001 HEIGHT  Height (cm)               64.7 WEEK 24 
#>  4 01-101-1001 HEIGHT  Height (cm)               68.2 WEEK 48 
#>  5 01-101-1001 WEIGHT  Weight (kg)                2.6 BASELINE
#>  6 01-101-1001 WEIGHT  Weight (kg)                5.3 WEEK 12 
#>  7 01-101-1001 WEIGHT  Weight (kg)                6.7 WEEK 24 
#>  8 01-101-1001 WEIGHT  Weight (kg)                7.4 WEEK 48 
#>  9 01-101-1001 BMI     Body Mass Index (kg/m^2)  11.7 BASELINE
#> 10 01-101-1001 BMI     Body Mass Index (kg/m^2)  15.2 WEEK 12 
#> 11 01-101-1001 BMI     Body Mass Index (kg/m^2)  16.0 WEEK 24 
#> 12 01-101-1001 BMI     Body Mass Index (kg/m^2)  15.9 WEEK 48 
```
