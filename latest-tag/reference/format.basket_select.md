# Returns a Character Representation of a `basket_select()` Object

The function returns a character representation of a
[`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
object. It can be used for error messages for example.

## Usage

``` r
# S3 method for class 'basket_select'
format(x, ...)
```

## Arguments

- x:

  A
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object

  Default value

  :   none

- ...:

  Not used

  Default value

  :   none

## Value

A character representation of the
[`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
object

## See also

[`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/admiral-package.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/extract_duplicate_records.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/signal_duplicate_records.md)

## Examples

``` r
format(basket_select(id = 42, scope = "NARROW", type = "smq"))
#> [1] "basket_select(name = NULL, id = 42, scope = \"NARROW\", type = \"smq\")"
```
