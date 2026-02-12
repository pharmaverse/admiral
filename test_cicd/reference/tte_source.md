# Create a `tte_source` Object

The `tte_source` object is used to define events and possible
censorings.

## Usage

``` r
tte_source(
  dataset_name,
  filter = NULL,
  date,
  censor = 0,
  set_values_to = NULL,
  order = order
)
```

## Arguments

- dataset_name:

  The name of the source dataset

  The name refers to the dataset provided by the `source_datasets`
  parameter of
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_tte.md).

  Default value

  :   none

- filter:

  An unquoted condition for selecting the observations from `dataset`
  which are events or possible censoring time points.

  Default value

  :   `NULL`

- date:

  A variable or expression providing the date of the event or censoring.
  A date, or a datetime can be specified. An unquoted symbol or
  expression is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  or
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/convert_dtc_to_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   none

- censor:

  Censoring value

  CDISC strongly recommends using `0` for events and positive integers
  for censoring.

  Default value

  :   `0`

- set_values_to:

  A named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining the variables to be set for the event or censoring, e.g.
  `exprs(EVENTDESC = "DEATH", SRCDOM = "ADSL", SRCVAR = "DTHDT")`. The
  values must be a symbol, a character string, a numeric value, an
  expression, or `NA`.

  Default value

  :   `NULL`

- order:

  Sort order

  An optional named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining additional variables that the source dataset is sorted on
  after `date`.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
      e.g. `exprs(ASEQ)`.

  Default value

  :   `order`

## Value

An object of class `tte_source`

## See also

[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_tte.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/test_cicd/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/records_source.md)
