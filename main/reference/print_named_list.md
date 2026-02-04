# Print Named List

Print Named List

## Usage

``` r
print_named_list(list, indent = 0)
```

## Arguments

- list:

  A named list

  Default value

  :   none

- indent:

  Indent

  The output is indented by the specified number of characters.

  Default value

  :   `0`

## Value

No return value, called for side effects

## See also

Utilities for printing:
[`print.adam_templates()`](https:/pharmaverse.github.io/admiral/main/reference/print.adam_templates.md),
[`print.duplicates()`](https:/pharmaverse.github.io/admiral/main/reference/print.duplicates.md),
[`print.source()`](https:/pharmaverse.github.io/admiral/main/reference/print.source.md)

## Examples

``` r
print_named_list(death_event)
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
