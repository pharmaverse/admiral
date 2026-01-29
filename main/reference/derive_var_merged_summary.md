# Merge Summary Variables

**\[deprecated\]** The `derive_var_merged_summary()` function has been
deprecated in favor of
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md).

## Usage

``` r
derive_var_merged_summary(
  dataset,
  dataset_add,
  by_vars,
  new_vars = NULL,
  filter_add = NULL,
  missing_values = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_add:

  Additional dataset

  The variables specified by the `by_vars` and the variables used on the
  left hand sides of the `new_vars` arguments are expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  The expressions on the left hand sides of `new_vars` are evaluated by
  the specified *variables*. Then the resulting values are merged to the
  input dataset (`dataset`) by the specified *variables*.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- new_vars:

  New variables to add

  The specified variables are added to the input dataset.

  A named list of expressions is expected:

  - LHS refer to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, an expression or NA. If summary
    functions are used, the values are summarized by the variables
    specified for `by_vars`. Any expression on the RHS must result in a
    single value per by group.

  For example:

        new_vars = exprs(
          DOSESUM = sum(AVAL),
          DOSEMEAN = mean(AVAL)
        )

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations fulfilling the specified condition are taken into
  account for summarizing. If the argument is not specified, all
  observations are considered.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- missing_values:

  Values for non-matching observations

  For observations of the input dataset (`dataset`) which do not have a
  matching observation in the additional dataset (`dataset_add`) the
  values of the specified variables are set to the specified value. Only
  variables specified for `new_vars` can be specified for
  `missing_values`.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   `NULL`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variables specified for `new_vars`.

## Details

1.  The records from the additional dataset (`dataset_add`) are
    restricted to those matching the `filter_add` condition.

2.  The new variables (`new_vars`) are created for each by group
    (`by_vars`) in the additional dataset (`dataset_add`) by calling
    [`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html).
    I.e., all observations of a by group are summarized to a single
    observation.

3.  The new variables are merged to the input dataset. For observations
    without a matching observation in the additional dataset the new
    variables are set to `NA`. Observations in the additional dataset
    which have no matching observation in the input dataset are ignored.

## See also

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
