# Create a `censor_source` Object

`censor_source` objects are used to define censorings as input for the
[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
function.

**Note:** This is a wrapper function for the more generic
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md).

## Usage

``` r
censor_source(
  dataset_name,
  filter = NULL,
  date,
  censor = 1,
  set_values_to = NULL,
  order = NULL
)
```

## Arguments

- dataset_name:

  The name of the source dataset

  The name refers to the dataset provided by the `source_datasets`
  parameter of
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md).

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
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  or
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md)
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
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  defining the variables to be set for the event or censoring, e.g.
  `exprs(EVENTDESC = "DEATH", SRCDOM = "ADSL", SRCVAR = "DTHDT")`. The
  values must be a symbol, a character string, a numeric value, an
  expression, or `NA`.

  Default value

  :   `NULL`

- order:

  Sort order

  An optional named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  defining additional variables that the source dataset is sorted on
  after `date`.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
      e.g. `exprs(ASEQ)`.

  Default value

  :   `order`

## Value

An object of class `censor_source`, inheriting from class `tte_source`

## See also

[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md),
[`death_event`](https:/pharmaverse.github.io/admiral/main/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/main/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/main/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)

## Examples

``` r
# Last study date known alive censor

censor_source(
  dataset_name = "adsl",
  date = LSTALVDT,
  set_values_to = exprs(
    EVNTDESC = "ALIVE",
    SRCDOM = "ADSL",
    SRCVAR = "LSTALVDT"
  )
)
#> <censor_source> object
#> dataset_name: "adsl"
#> filter: NULL
#> date: LSTALVDT
#> censor: 1
#> set_values_to:
#>   EVNTDESC: "ALIVE"
#>   SRCDOM: "ADSL"
#>   SRCVAR: "LSTALVDT"
#> order: NULL
```
