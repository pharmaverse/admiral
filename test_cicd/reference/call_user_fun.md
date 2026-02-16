# Calls a Function Provided by the User

**\[deprecated\]**

Calls a function provided by the user and adds the function call to the
error message if the call fails.

## Usage

``` r
call_user_fun(call)
```

## Arguments

- call:

  Call to be executed

  Default value

  :   none

## Value

The return value of the function call

## See also

Other deprecated:
[`date_source()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/test_cicd/reference/get_summary_records.md)

## Examples

``` r
call_user_fun(compute_bmi(
  height = 172,
  weight = 60
))
#> [1] 20.28123

try(call_user_fun(compute_bmi(
  height = 172,
  weight = "hallo"
)))
#> Error in call_user_fun(compute_bmi(height = 172, weight = "hallo")) : 
#>   Calling `compute_bmi(height = 172, weight = "hallo")` caused the
#> following error:
#> Argument `weight` must be a numeric vector, but it is a string.
```
