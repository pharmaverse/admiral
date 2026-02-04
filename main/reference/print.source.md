# Print `source` Objects

Print `source` Objects

## Usage

``` r
# S3 method for class 'source'
print(x, ...)
```

## Arguments

- x:

  An `source` object

  Default value

  :   none

- ...:

  If `indent = <numeric value>` is specified the output is indented by
  the specified number of characters.

  Default value

  :   none

## Value

No return value, called for side effects

## See also

Utilities for printing:
[`print.adam_templates()`](https:/pharmaverse.github.io/admiral/main/reference/print.adam_templates.md),
[`print.duplicates()`](https:/pharmaverse.github.io/admiral/main/reference/print.duplicates.md),
[`print_named_list()`](https:/pharmaverse.github.io/admiral/main/reference/print_named_list.md)

## Examples

``` r
print(death_event)
#> <event_source> object
#> dataset_name: "adsl"
#> filter: DTHFL == "Y"
#> date: DTHDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "DEATH"
#>   SRCDOM: "ADSL"
#>   SRCVAR: "DTHDT"
#> order: NULL
```
