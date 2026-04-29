# Derive Expected Records

Add expected records as new observations for each 'by group' when the
dataset contains missing observations.

## Usage

``` r
derive_expected_records(
  dataset,
  dataset_ref,
  by_vars = NULL,
  set_values_to = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `dataset_ref` and `by_vars` arguments
  are expected to be in the dataset.

  Default value

  :   none

- dataset_ref:

  Expected observations dataset

  Data frame with the expected observations, e.g., all the expected
  combinations of `PARAMCD`, `PARAM`, `AVISIT`, `AVISITN`, ...

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` those observations from
  `dataset_ref` are added to the output dataset which do not have a
  corresponding observation in the input dataset.

  Default value

  :   `NULL`

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations.

  A list of variable name-value pairs is expected.

  - LHS refers to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, `NA`, or expressions, e.g.,
    `exprs(PARAMCD = "TDOSE", PARCAT1 = "OVERALL")`.

  Default value

  :   `NULL`

## Value

The input dataset with the missed expected observations added for each
`by_vars`. Note, a variable will only be populated in the new parameter
rows if it is specified in `by_vars` or `set_values_to`.

## Details

For each group (the variables specified in the `by_vars` parameter),
those records from `dataset_ref` that are missing in the input dataset
are added to the output dataset.

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/default_qtc_paramcd.md),
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
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

adqs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
  "1",      "a",             1, "WEEK 1",   10,
  "1",      "b",             1, "WEEK 1",   11,
  "2",      "a",             2, "WEEK 2",   12,
  "2",      "b",             2, "WEEK 2",   14
)

# Example 1. visit variables are parameter independent
parm_visit_ref <- tribble(
  ~AVISITN, ~AVISIT,
  1,        "WEEK 1",
  2,        "WEEK 2"
)

derive_expected_records(
  dataset = adqs,
  dataset_ref = parm_visit_ref,
  by_vars = exprs(USUBJID, PARAMCD),
  set_values_to = exprs(DTYPE = "DERIVED")
)
#> # A tibble: 8 × 6
#>   USUBJID PARAMCD AVISITN AVISIT  AVAL DTYPE  
#>   <chr>   <chr>     <dbl> <chr>  <dbl> <chr>  
#> 1 1       a             1 WEEK 1    10 NA     
#> 2 1       a             2 WEEK 2    NA DERIVED
#> 3 1       b             1 WEEK 1    11 NA     
#> 4 1       b             2 WEEK 2    NA DERIVED
#> 5 2       a             1 WEEK 1    NA DERIVED
#> 6 2       a             2 WEEK 2    12 NA     
#> 7 2       b             1 WEEK 1    NA DERIVED
#> 8 2       b             2 WEEK 2    14 NA     

# Example 2. visit variables are parameter dependent
parm_visit_ref <- tribble(
  ~PARAMCD, ~AVISITN, ~AVISIT,
  "a",             1, "WEEK 1",
  "a",             2, "WEEK 2",
  "b",             1, "WEEK 1"
)

derive_expected_records(
  dataset = adqs,
  dataset_ref = parm_visit_ref,
  by_vars = exprs(USUBJID, PARAMCD),
  set_values_to = exprs(DTYPE = "DERIVED")
)
#> # A tibble: 7 × 6
#>   USUBJID PARAMCD AVISITN AVISIT  AVAL DTYPE  
#>   <chr>   <chr>     <dbl> <chr>  <dbl> <chr>  
#> 1 1       a             1 WEEK 1    10 NA     
#> 2 1       a             2 WEEK 2    NA DERIVED
#> 3 1       b             1 WEEK 1    11 NA     
#> 4 2       a             1 WEEK 1    NA DERIVED
#> 5 2       a             2 WEEK 2    12 NA     
#> 6 2       b             1 WEEK 1    NA DERIVED
#> 7 2       b             2 WEEK 2    14 NA     
```
