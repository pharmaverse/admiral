# Create a `basket_select` object

Create a `basket_select` object

## Usage

``` r
basket_select(name = NULL, id = NULL, scope = NULL, type, ...)
```

## Arguments

- name:

  Name of the query used to select the definition of the query from the
  company database.

  Default value

  :   `NULL`

- id:

  Identifier of the query used to select the definition of the query
  from the company database.

  Default value

  :   `NULL`

- scope:

  Scope of the query used to select the definition of the query from the
  company database.

  Permitted values

  :   `"BROAD"`, `"NARROW"`, `NA_character_`

  Default value

  :   `NULL`

- type:

  The type argument expects a character scalar. It is passed to the
  company specific get_terms() function such that the function can
  determine which sort of basket is requested

  Default value

  :   none

- ...:

  Any number of *named* function arguments. Can be used to pass in
  company specific conditions or flags that will then be used in
  user-defined function that is passed into argument `get_terms_fun` for
  function
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md).

  Default value

  :   none

## Value

An object of class `basket_select`.

## Details

Exactly one of `name` or `id` must be specified.

## See also

[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md),
[`query()`](https:/pharmaverse.github.io/admiral/main/reference/query.md)

Source Objects:
[`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/main/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/main/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/main/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)
