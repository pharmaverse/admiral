# Map `"Y"` and `"N"` to Numeric Values

Map `"Y"` and `"N"` to numeric values.

## Usage

``` r
yn_to_numeric(arg)
```

## Arguments

- arg:

  Character vector

  Default value

  :   none

## Value

`1` if `arg` equals `"Y"`, `0` if `arg` equals `"N"`, `NA_real_`
otherwise

## See also

Utilities for Formatting Observations:
[`convert_blanks_to_na()`](https:/pharmaverse.github.io/admiral/main/reference/convert_blanks_to_na.md),
[`convert_na_to_blanks()`](https:/pharmaverse.github.io/admiral/main/reference/convert_na_to_blanks.md)

## Examples

``` r
yn_to_numeric(c("Y", "N", NA_character_))
#> [1]  1  0 NA
```
