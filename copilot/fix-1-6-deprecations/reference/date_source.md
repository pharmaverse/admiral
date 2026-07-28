# Create a `date_source` object

**\[deprecated\]** The `date_source()` function has been deprecated in
favor of
[`event()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/event.md).

Create a `date_source` object as input for
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dt.md)
and
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dtm.md).

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

[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dtm.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dt.md)

Other deprecated: `call_user_fun()`,
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/copilot/fix-1-6-deprecations/reference/get_summary_records.md)
