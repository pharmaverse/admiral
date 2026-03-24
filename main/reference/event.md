# Create a `event` Object

The `event` object is used to define events as input for the
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
and
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)
functions.

## Usage

``` r
event(
  dataset_name = NULL,
  condition = NULL,
  mode = NULL,
  order = NULL,
  set_values_to = NULL,
  keep_source_vars = NULL,
  description = NULL
)
```

## Arguments

- dataset_name:

  Dataset name of the dataset to be used as input for the event. The
  name refers to the dataset specified for `source_datasets` in
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md).
  If the argument is not specified, the input dataset (`dataset`) of
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  is used.

  Permitted values

  :   a character scalar

  Default value

  :   `NULL`

- condition:

  An unquoted condition for selecting the observations, which will
  contribute to the extreme event. If the condition contains summary
  functions like [`all()`](https://rdrr.io/r/base/all.html), they are
  evaluated for each by group separately.

  Permitted values

  :   an unquoted condition

  Default value

  :   `NULL`

- mode:

  If specified, the first or last observation with respect to `order` is
  selected for each by group.

  Permitted values

  :   `"first"`, `"last"`, `NULL`

  Default value

  :   `NULL`

- order:

  The specified variables or expressions are used to select the first or
  last observation if `mode` is specified.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))` or `NULL`

  Default value

  :   `NULL`

- set_values_to:

  A named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  defining the variables to be set for the event, e.g.
  `exprs(PARAMCD = "WSP", PARAM = "Worst Sleeping Problems")`. The
  values can be a symbol, a character string, a numeric value, `NA` or
  an expression.

  Permitted values

  :   a named list of expressions, e.g., created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)

  Default value

  :   `NULL`

- keep_source_vars:

  Variables to keep from the source dataset

  The specified variables are kept for the selected observations. The
  variables specified for `by_vars` (of
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md))
  and created by `set_values_to` are always kept.

  Permitted values

  :   A list of expressions where each element is a symbol or a
      tidyselect expression, e.g.,
      `exprs(VISIT, VISITNUM, starts_with("RS"))`.

  Default value

  :   `NULL`

- description:

  Description of the event

  The description does not affect the derivations where the event is
  used. It is intended for documentation only.

  Permitted values

  :   a character scalar

  Default value

  :   `NULL`

## Value

An object of class `event`

## See also

[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/main/reference/tte_source_objects.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/main/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/main/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)
