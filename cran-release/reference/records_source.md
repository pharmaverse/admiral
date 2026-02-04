# Create a `records_source` Object

The `records_source` object is used to find extreme records of interest.

## Usage

``` r
records_source(dataset_name, filter = NULL, new_vars)
```

## Arguments

- dataset_name:

  The name of the source dataset

  The name refers to the dataset provided by the `source_datasets`
  argument of
  [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_extreme_record.md).

  Default value

  :   none

- filter:

  An unquoted condition for selecting the observations from `dataset`.

  Default value

  :   `NULL`

- new_vars:

  Variables to add

  The specified variables from the source datasets are added to the
  output dataset. Variables can be renamed by naming the element, i.e.,
  `new_vars = exprs(<new name> = <old name>)`.

  For example `new_vars = exprs(var1, var2)` adds variables `var1` and
  `var2` from to the input dataset.

  And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
  `old_var2` from the source dataset and adds them to the input dataset
  renaming `old_var2` to `new_var2`. Expressions can be used to create
  new variables (see for example `new_vars` argument in
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))`

  Default value

  :   none

## Value

An object of class `records_source`

## See also

[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_extreme_record.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/query.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/tte_source.md)
