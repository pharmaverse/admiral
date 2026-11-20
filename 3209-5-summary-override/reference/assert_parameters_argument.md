# Asserts `parameters` Argument and Converts to List of Expressions

The function asserts that the argument is a character vector or a list
of expressions. If it is a character vector, it converts it to a list of
symbols.

## Usage

``` r
assert_parameters_argument(parameters, optional = TRUE)
```

## Arguments

- parameters:

  The argument to check

  Default value

  :   none

- optional:

  Is the checked argument optional? If set to `FALSE` and `parameters`
  is `NULL` then an error is thrown.

  Default value

  :   `TRUE`

## Value

The `parameters` argument (converted to a list of symbol, if it is a
character vector)
