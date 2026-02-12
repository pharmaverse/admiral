# Print `adam_templates` Objects

Print `adam_templates` Objects

## Usage

``` r
# S3 method for class 'adam_templates'
print(x, ...)
```

## Arguments

- x:

  A `adam_templates` object

  Default value

  :   none

- ...:

  Not used

  Default value

  :   none

## Value

No return value, called for side effects

## See also

[`list_all_templates()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/list_all_templates.md)

Utilities for printing:
[`print.duplicates()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/print.duplicates.md),
[`print.source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/print.source.md),
[`print_named_list()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/print_named_list.md)

## Examples

``` r
templates <- list_all_templates()
print(templates)
#> Existing ADaM templates in package 'admiral':
#> • ADAB
#> • ADAE
#> • ADCM
#> • ADEG
#> • ADEX
#> • ADLB
#> • ADLBHY
#> • ADMH
#> • ADPC
#> • ADPP
#> • ADPPK
#> • ADSL
#> • ADVS
```
