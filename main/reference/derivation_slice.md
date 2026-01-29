# Create a `derivation_slice` Object

Create a `derivation_slice` object as input for
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md).

## Usage

``` r
derivation_slice(filter, args = NULL)
```

## Arguments

- filter:

  An unquoted condition for defining the observations of the slice

  Default value

  :   none

- args:

  Arguments of the derivation to be used for the slice

  A
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
  object is expected.

  Default value

  :   `NULL`

## Value

An object of class `derivation_slice`

## See also

[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md),
[`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)

Higher Order Functions:
[`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md),
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)
