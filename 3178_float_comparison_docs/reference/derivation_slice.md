# Create a `derivation_slice` Object

Create a `derivation_slice` object as input for
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/slice_derivation.md).

## Usage

``` r
derivation_slice(filter, args = NULL)
```

## Arguments

- filter:

  An unquoted condition for defining the observations of the slice

  Comparing derived numeric variables to fixed values, e.g.,
  `PCHG <= -90`, may give unexpected results due to floating point
  representation. For details and solutions see the "Floating Point
  Comparisons" section in
  [`vignette("concepts_conventions")`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/articles/concepts_conventions.md).

  Default value

  :   none

- args:

  Arguments of the derivation to be used for the slice

  A
  [`params()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/params.md)
  object is expected.

  Default value

  :   `NULL`

## Value

An object of class `derivation_slice`

## See also

[`slice_derivation()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/slice_derivation.md),
[`params()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/params.md)

Higher Order Functions:
[`call_derivation()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/call_derivation.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/restrict_derivation.md),
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/3178_float_comparison_docs/reference/slice_derivation.md)
