# Adds a Parameter for Dose Intensity

Adds a record for the dose intensity for each by group (e.g., subject
and visit) where the source parameters are available.

**Note:** This is a wrapper function for the more generic
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md).

The analysis value of the new parameter is derived as Total Dose /
Planned Dose \* 100

## Usage

``` r
derive_param_doseint(
  dataset,
  by_vars,
  set_values_to = exprs(PARAMCD = "TNDOSINT"),
  tadm_code = "TNDOSE",
  tpadm_code = "TSNDOSE",
  zero_doses = "Inf",
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
  (`filter` parameter) and to the parameters specified by `tadm_code`
  and `padm_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  Only variables specified in `by_vars` will be populated in the newly
  created records.

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

- tadm_code:

  Total Doses Administered parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the total dose administered. The `AVAL` associated with
  this `PARAMCD` will be the numerator of the dose intensity
  calculation.

  Permitted values

  :   character value

  Default value

  :   `"TNDOSE"`

- tpadm_code:

  Total Doses Planned parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the total planned dose. The `AVAL` associated with this
  `PARAMCD` will be the denominator of the dose intensity calculation.

  Permitted values

  :   character value

  Default value

  :   `"TSNDOSE"`

- zero_doses:

  Flag indicating logic for handling 0 planned or administered doses for
  a `by_vars` group

  Permitted values

  :   `Inf`, `100`

      No record is returned if either the planned (`tpadm_code`) or
      administered (`tadm_code`) `AVAL` are `NA`. No record is returned
      is a record does not exist for both `tadm_code` and `tpadm_code`
      for the specified `by_var`.

      If `zero_doses` = `Inf`:

      1.  If the planned dose (`tpadm_code`) is 0 and administered dose
          (`tadm_code`) is 0, `NaN` is returned.

      2.  If the planned dose (`tpadm_code`) is 0 and the administered
          dose (`tadm_code`) is \> 0, `Inf` is returned.

      If `zero_doses` = `100` :

      1.  If the planned dose (`tpadm_code`) is 0 and administered dose
          (`tadm_code`) is 0, 0 is returned.

      2.  If the planned dose (`tpadm_code`) is 0 and the administered
          dose (`tadm_code`) is \> 0, 100 is returned.

  Default value

  :   `"Inf"`

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

The input dataset with the new parameter rows added. Note, a variable
will only be populated in the new parameter rows if it is specified in
`by_vars`.

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
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
library(lubridate, warn.conflicts = FALSE)

adex <- tribble(
  ~USUBJID, ~PARAMCD, ~VISIT, ~ANL01FL, ~ASTDT, ~AENDT, ~AVAL,
  "P001", "TNDOSE", "V1", "Y", ymd("2020-01-01"), ymd("2020-01-30"), 59,
  "P001", "TSNDOSE", "V1", "Y", ymd("2020-01-01"), ymd("2020-02-01"), 96,
  "P001", "TNDOSE", "V2", "Y", ymd("2020-02-01"), ymd("2020-03-15"), 88,
  "P001", "TSNDOSE", "V2", "Y", ymd("2020-02-05"), ymd("2020-03-01"), 88,
  "P002", "TNDOSE", "V1", "Y", ymd("2021-01-01"), ymd("2021-01-30"), 0,
  "P002", "TSNDOSE", "V1", "Y", ymd("2021-01-01"), ymd("2021-02-01"), 0,
  "P002", "TNDOSE", "V2", "Y", ymd("2021-02-01"), ymd("2021-03-15"), 52,
  "P002", "TSNDOSE", "V2", "Y", ymd("2021-02-05"), ymd("2021-03-01"), 0
)

derive_param_doseint(
  adex,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(PARAMCD = "TNDOSINT"),
  tadm_code = "TNDOSE",
  tpadm_code = "TSNDOSE"
)
#> # A tibble: 12 × 7
#>    USUBJID PARAMCD  VISIT ANL01FL ASTDT      AENDT       AVAL
#>    <chr>   <chr>    <chr> <chr>   <date>     <date>     <dbl>
#>  1 P001    TNDOSE   V1    Y       2020-01-01 2020-01-30  59  
#>  2 P001    TSNDOSE  V1    Y       2020-01-01 2020-02-01  96  
#>  3 P001    TNDOSE   V2    Y       2020-02-01 2020-03-15  88  
#>  4 P001    TSNDOSE  V2    Y       2020-02-05 2020-03-01  88  
#>  5 P002    TNDOSE   V1    Y       2021-01-01 2021-01-30   0  
#>  6 P002    TSNDOSE  V1    Y       2021-01-01 2021-02-01   0  
#>  7 P002    TNDOSE   V2    Y       2021-02-01 2021-03-15  52  
#>  8 P002    TSNDOSE  V2    Y       2021-02-05 2021-03-01   0  
#>  9 P001    TNDOSINT V1    NA      NA         NA          61.5
#> 10 P001    TNDOSINT V2    NA      NA         NA         100  
#> 11 P002    TNDOSINT V1    NA      NA         NA         NaN  
#> 12 P002    TNDOSINT V2    NA      NA         NA         Inf  

derive_param_doseint(
  adex,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(PARAMCD = "TDOSINT2"),
  tadm_code = "TNDOSE",
  tpadm_code = "TSNDOSE",
  zero_doses = "100"
)
#> # A tibble: 12 × 7
#>    USUBJID PARAMCD  VISIT ANL01FL ASTDT      AENDT       AVAL
#>    <chr>   <chr>    <chr> <chr>   <date>     <date>     <dbl>
#>  1 P001    TNDOSE   V1    Y       2020-01-01 2020-01-30  59  
#>  2 P001    TSNDOSE  V1    Y       2020-01-01 2020-02-01  96  
#>  3 P001    TNDOSE   V2    Y       2020-02-01 2020-03-15  88  
#>  4 P001    TSNDOSE  V2    Y       2020-02-05 2020-03-01  88  
#>  5 P002    TNDOSE   V1    Y       2021-01-01 2021-01-30   0  
#>  6 P002    TSNDOSE  V1    Y       2021-01-01 2021-02-01   0  
#>  7 P002    TNDOSE   V2    Y       2021-02-01 2021-03-15  52  
#>  8 P002    TSNDOSE  V2    Y       2021-02-05 2021-03-01   0  
#>  9 P001    TDOSINT2 V1    NA      NA         NA          61.5
#> 10 P001    TDOSINT2 V2    NA      NA         NA         100  
#> 11 P002    TDOSINT2 V1    NA      NA         NA           0  
#> 12 P002    TDOSINT2 V2    NA      NA         NA         100  
```
