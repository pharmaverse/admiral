# Add a parameter for lab differentials converted to absolute values

Add a parameter by converting lab differentials from fraction or
percentage to absolute values

## Usage

``` r
derive_param_wbc_abs(
  dataset,
  by_vars,
  set_values_to,
  get_unit_expr,
  wbc_unit = "10^9/L",
  wbc_code = "WBC",
  diff_code,
  diff_type = "fraction"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD`, and `AVAL` are expected as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset, and to the parameters specified by `wbc_code`
  and `diff_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   none

- set_values_to:

  Variables to set

  A named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining the variables to be set for the new parameter, e.g.
  `exprs(PARAMCD = "LYMPH", PARAM = "Lymphocytes Abs (10^9/L)")` is
  expected.

  Default value

  :   none

- get_unit_expr:

  An expression providing the unit of the parameter

  The result is used to check the units of the input parameters.

  Permitted values

  :   a variable containing unit from the input dataset, or a function
      call, for example, `get_unit_expr = extract_unit(PARAM)`.

  Default value

  :   none

- wbc_unit:

  A string containing the required unit of the WBC parameter

  Default value

  :   `"10^9/L"`

- wbc_code:

  White Blood Cell (WBC) parameter

  The observations where `PARAMCD` equals the specified value are
  considered as the WBC absolute results to use for converting the
  differentials.

  Permitted values

  :   character value

  Default value

  :   `"WBC"`

- diff_code:

  white blood differential parameter

  The observations where `PARAMCD` equals the specified value are
  considered as the white blood differential lab results in fraction or
  percentage value to be converted into absolute value.

  Default value

  :   none

- diff_type:

  A string specifying the type of differential

  Permitted values

  :   `"percent"`, `"fraction"`

  Default value

  :   `"fraction"`

## Value

The input dataset with the new parameter added

## Details

If `diff_type` is `"percent"`, the analysis value of the new parameter
is derived as \$\$\frac{White Blood Cell Count \* Percentage
Value}{100}\$\$

If `diff_type` is `"fraction"`, the analysis value of the new parameter
is derived as \$\$White Blood Cell Count \* Fraction Value\$\$

New records are created for each group of records (grouped by `by_vars`)
if 1) the white blood cell component in absolute value is not already
available from the input dataset, and 2) the white blood cell absolute
value (identified by `wbc_code`) and the white blood cell differential
(identified by `diff_code`) are both present.

## See also

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
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_rr.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

test_lb <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
  "P01", "WBC", 33, "Leukocyte Count (10^9/L)", "CYCLE 1 DAY 1",
  "P01", "WBC", 38, "Leukocyte Count (10^9/L)", "CYCLE 2 DAY 1",
  "P01", "LYMLE", 0.90, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
  "P01", "LYMLE", 0.70, "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
  "P01", "ALB", 36, "Albumin (g/dL)", "CYCLE 2 DAY 1",
  "P02", "WBC", 33, "Leukocyte Count (10^9/L)", "CYCLE 1 DAY 1",
  "P02", "LYMPH", 29, "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1",
  "P02", "LYMLE", 0.87, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
  "P03", "LYMLE", 0.89, "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1"
)

derive_param_wbc_abs(
  dataset = test_lb,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(
    PARAMCD = "LYMPH",
    PARAM = "Lymphocytes Abs (10^9/L)",
    DTYPE = "CALCULATION"
  ),
  get_unit_expr = extract_unit(PARAM),
  wbc_code = "WBC",
  diff_code = "LYMLE",
  diff_type = "fraction"
)
#> # A tibble: 11 × 6
#>    USUBJID PARAMCD  AVAL PARAM                       VISIT         DTYPE      
#>    <chr>   <chr>   <dbl> <chr>                       <chr>         <chr>      
#>  1 P01     WBC     33    Leukocyte Count (10^9/L)    CYCLE 1 DAY 1 NA         
#>  2 P01     WBC     38    Leukocyte Count (10^9/L)    CYCLE 2 DAY 1 NA         
#>  3 P01     LYMLE    0.9  Lymphocytes (fraction of 1) CYCLE 1 DAY 1 NA         
#>  4 P01     LYMLE    0.7  Lymphocytes (fraction of 1) CYCLE 2 DAY 1 NA         
#>  5 P01     ALB     36    Albumin (g/dL)              CYCLE 2 DAY 1 NA         
#>  6 P02     WBC     33    Leukocyte Count (10^9/L)    CYCLE 1 DAY 1 NA         
#>  7 P02     LYMPH   29    Lymphocytes Abs (10^9/L)    CYCLE 1 DAY 1 NA         
#>  8 P02     LYMLE    0.87 Lymphocytes (fraction of 1) CYCLE 1 DAY 1 NA         
#>  9 P03     LYMLE    0.89 Lymphocytes (fraction of 1) CYCLE 1 DAY 1 NA         
#> 10 P01     LYMPH   29.7  Lymphocytes Abs (10^9/L)    CYCLE 1 DAY 1 CALCULATION
#> 11 P01     LYMPH   26.6  Lymphocytes Abs (10^9/L)    CYCLE 2 DAY 1 CALCULATION
```
