# Create a `derivation_slice` Object

Create a `derivation_slice` object as input for
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/slice_derivation.md).

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
  [`params()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/params.md)
  object is expected.

  Default value

  :   `NULL`

## Value

An object of class `derivation_slice`

## See also

[`slice_derivation()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/slice_derivation.md),
[`params()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/params.md)

Higher Order Functions:
[`call_derivation()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/call_derivation.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/restrict_derivation.md),
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/cran-release/3006_template_ids_and_links/reference/slice_derivation.md)
