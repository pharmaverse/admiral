# Turn a Character Vector into a List of Expressions

Turn a character vector into a list of expressions

## Usage

``` r
chr2vars(chr)
```

## Arguments

- chr:

  A character vector

  Default value

  :   none

## Value

A `list` of expressions as returned by
[`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md)

## See also

Utilities for working with quosures/list of expressions:
[`negate_vars()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/negate_vars.md)

## Examples

``` r
chr2vars(c("USUBJID", "AVAL"))
#> [[1]]
#> USUBJID
#> 
#> [[2]]
#> AVAL
#> 
```
