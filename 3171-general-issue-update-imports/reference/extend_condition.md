# Extend a condition string by adding a new condition based on a variable and its value

This internal helper function extends a condition string by appending a
new condition that checks if a variable equals a specific value.

## Usage

``` r
extend_condition(cond, var, is)
```

## Arguments

- cond:

  A character string representing an existing condition.

  Default value

  :   none

- var:

  A character string representing the name of the variable to check.

  Default value

  :   none

- is:

  A character string representing the value the variable should be equal
  to.

  Default value

  :   none

## Value

A character string representing the extended condition.

## Examples

``` r
# Extend an existing condition to include a check for 'AGE == "30"'
admiral:::extend_condition("SEX == 'M'", "AGE", "30")
#> [1] "SEX == 'M' & AGE == '30'"
```
