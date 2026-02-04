# Create a `dthcaus_source` Object

**\[deprecated\]** The `dthcaus_source()` function and
`dthcaus_source()` have been deprecated in favor of
[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md).

## Usage

``` r
dthcaus_source(
  dataset_name,
  filter,
  date,
  order = NULL,
  mode = "first",
  dthcaus,
  set_values_to = NULL
)
```

## Arguments

- dataset_name:

  The name of the dataset, i.e. a string, used to search for the death
  cause.

  Default value

  :   none

- filter:

  An expression used for filtering `dataset`.

  Default value

  :   none

- date:

  A date or datetime variable or an expression to be used for sorting
  `dataset`.

  Default value

  :   none

- order:

  Sort order

  Additional variables/expressions to be used for sorting the `dataset`.
  The dataset is ordered by `date` and `order`. Can be used to avoid
  duplicate record warning.

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))` or `NULL`

  Default value

  :   `NULL`

- mode:

  One of `"first"` or `"last"`. Either the `"first"` or `"last"`
  observation is preserved from the `dataset` which is ordered by
  `date`.

  Default value

  :   `"first"`

- dthcaus:

  A variable name, an expression, or a string literal

  If a variable name is specified, e.g., `AEDECOD`, it is the variable
  in the source dataset to be used to assign values to `DTHCAUS`; if an
  expression, e.g., `str_to_upper(AEDECOD)`, it is evaluated in the
  source dataset and the results is assigned to `DTHCAUS`; if a string
  literal, e.g. `"Adverse Event"`, it is the fixed value to be assigned
  to `DTHCAUS`.

  Default value

  :   none

- set_values_to:

  Variables to be set to trace the source dataset

  Default value

  :   `NULL`

## Value

An object of class "dthcaus_source".

## See also

[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_dthcaus.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_summary.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_summary_records.md)

## Examples

``` r
# Deaths sourced from AE
src_ae <- dthcaus_source(
  dataset_name = "ae",
  filter = AEOUT == "FATAL",
  date = AEDTHDT,
  mode = "first",
  dthcaus = AEDECOD
)

# Deaths sourced from DS
src_ds <- dthcaus_source(
  dataset_name = "ds",
  filter = DSDECOD == "DEATH",
  date = convert_dtc_to_dt(DSSTDTC),
  mode = "first",
  dthcaus = DSTERM
)
```
