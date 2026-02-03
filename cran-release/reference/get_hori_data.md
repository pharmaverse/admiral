# Creating Temporary Parameters and `<variable>.<parameter>` Variables

The function creates temporary parameters and variables of the form
`<variable>.<parameter>`, e.g., `AVAL.WEIGHT`.

## Usage

``` r
get_hori_data(dataset, by_vars, parameters, set_values_to, filter)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   none

- parameters:

  List of parameter codes

  The input dataset is restricted to the specified parameter codes. If
  an expression is specified, a new parameter code is added to the input
  dataset. The name of the element defines the parameter code and the
  expression the observations to select.

  Permitted values

  :   A character vector of `PARAMCD` values or a list of expressions

  Default value

  :   none

- set_values_to:

  All variables of the form `<variable>.<parameter>` like `AVAL.WEIGHT`
  are added to the input dataset. They are set to the value of the
  variable for the parameter. E.g., `AVAL.WEIGHT` is set to the value of
  `AVAL` where `PARAMCD == "WEIGHT"`.

  Permitted values

  :   

  Default value

  :   none

- filter:

  Filter condition used for restricting the input dataset

  The specified filter condition is used in the warnings only. It is not
  applied to the input dataset.

  Permitted values

  :   An unquoted expression

  Default value

  :   none

## Value

A dataset with one observation per by group. It contains the variables
specified for `by_vars` and all variables of the form
`<variable>.<parameter>` occurring in `set_values_to`.
