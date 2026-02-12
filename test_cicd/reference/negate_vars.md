# Negate List of Variables

The function adds a minus sign as prefix to each variable.

## Usage

``` r
negate_vars(vars = NULL)
```

## Arguments

- vars:

  List of variables created by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)

  Default value

  :   `NULL`

## Value

A list of expressions

## Details

This is useful if a list of variables should be removed from a dataset,
e.g., `select(!!!negate_vars(by_vars))` removes all by variables.

## See also

Utilities for working with quosures/list of expressions:
[`chr2vars()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/chr2vars.md)

## Examples

``` r
negate_vars(exprs(USUBJID, STUDYID))
#> [[1]]
#> -USUBJID
#> 
#> [[2]]
#> -STUDYID
#> 
```
