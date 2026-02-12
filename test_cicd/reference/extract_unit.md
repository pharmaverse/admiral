# Extract Unit From Parameter Description

Extract the unit of a parameter from a description like "Param (unit)".

## Usage

``` r
extract_unit(x)
```

## Arguments

- x:

  A parameter description

  Default value

  :   none

## Value

A string

## See also

Utilities used within Derivation functions:
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_vars_query.md)

## Examples

``` r
extract_unit("Height (cm)")
#> [1] "cm"

extract_unit("Diastolic Blood Pressure (mmHg)")
#> [1] "mmHg"
```
