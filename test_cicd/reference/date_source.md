# Create a `date_source` object

**\[deprecated\]** The `date_source()` function has been deprecated in
favor of
[`event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event.md).

Create a `date_source` object as input for
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dt.md)
and
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md).

## Usage

``` r
date_source(dataset_name, filter = NULL, date, set_values_to = NULL)
```

## Arguments

- dataset_name:

  The name of the dataset, i.e. a string, used to search for the date.

  Default value

  :   none

- filter:

  An unquoted condition for filtering `dataset`.

  Default value

  :   `NULL`

- date:

  A variable or an expression providing a date. A date or a datetime can
  be specified. An unquoted symbol or expression is expected.

  Default value

  :   none

- set_values_to:

  Variables to be set

  Default value

  :   `NULL`

## Value

An object of class `date_source`.

## See also

[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dt.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/call_user_fun.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_summary_records.md)

## Examples

``` r
# treatment end date from ADSL
trt_end_date <- date_source(
  dataset_name = "adsl",
  date = TRTEDT
)
#> Warning: `date_source()` was deprecated in admiral 1.2.0.
#> ℹ Please use `event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation

# lab date from LB where assessment was taken, i.e. not "NOT DONE"
lb_date <- date_source(
  dataset_name = "lb",
  filter = LBSTAT != "NOT DONE" | is.na(LBSTAT),
  date = convert_dtc_to_dt(LBDTC)
)

# death date from ADSL including traceability variables
death_date <- date_source(
  dataset_name = "adsl",
  date = DTHDT,
  set_values_to = exprs(
    LALVDOM = "ADSL",
    LALVVAR = "DTHDT"
  )
)
```
