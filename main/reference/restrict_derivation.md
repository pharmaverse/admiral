# Execute a Derivation on a Subset of the Input Dataset

Execute a derivation on a subset of the input dataset.

## Usage

``` r
restrict_derivation(dataset, derivation, args = NULL, filter)
```

## Arguments

- dataset:

  Input dataset

  Default value

  :   none

- derivation:

  Derivation

  A function that performs a specific derivation is expected. A
  derivation adds variables or observations to a dataset. The first
  argument of a derivation must expect a dataset and the derivation must
  return a dataset. All expected arguments for the derivation function
  must be provided through the
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
  objects passed to the `args` argument.

  Default value

  :   none

- args:

  Arguments of the derivation

  A
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
  object is expected.

  Default value

  :   `NULL`

- filter:

  Filter condition

  Default value

  :   none

## Details

It is also possible to pass functions from outside the `{admiral}`
package to `restrict_derivation()`, e.g. an extension package function,
or
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
The only requirement for a function being passed to `derivation` is that
it must take a dataset as its first argument and return a dataset.

## See also

[`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)
[`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md)

Higher Order Functions:
[`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md),
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md),
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)

## Examples

``` r
library(tibble)

adlb <- tribble(
  ~USUBJID, ~AVISITN, ~AVAL, ~ABLFL,
  "1",            -1,   113, NA_character_,
  "1",             0,   113, "Y",
  "1",             3,   117, NA_character_,
  "2",             0,    95, "Y",
  "3",             0,   111, "Y",
  "3",             1,   101, NA_character_,
  "3",             2,   123, NA_character_
)

# Derive BASE for post-baseline records only (derive_var_base() can not be used in this case
# as it requires the baseline observation to be in the input dataset)
restrict_derivation(
  adlb,
  derivation = derive_vars_merged,
  args = params(
    by_vars = exprs(USUBJID),
    dataset_add = adlb,
    filter_add = ABLFL == "Y",
    new_vars = exprs(BASE = AVAL)
  ),
  filter = AVISITN > 0
)
#> # A tibble: 7 × 5
#>   USUBJID AVISITN  AVAL ABLFL  BASE
#>   <chr>     <dbl> <dbl> <chr> <dbl>
#> 1 1             3   117 NA      113
#> 2 3             1   101 NA      111
#> 3 3             2   123 NA      111
#> 4 1            -1   113 NA       NA
#> 5 1             0   113 Y        NA
#> 6 2             0    95 Y        NA
#> 7 3             0   111 Y        NA

# Derive BASE for baseline and post-baseline records only
restrict_derivation(
  adlb,
  derivation = derive_var_base,
  args = params(
    by_vars = exprs(USUBJID)
  ),
  filter = AVISITN >= 0
) %>%
  # Derive CHG for post-baseline records only
  restrict_derivation(
    derivation = derive_var_chg,
    filter = AVISITN > 0
  )
#> # A tibble: 7 × 6
#>   USUBJID AVISITN  AVAL ABLFL  BASE   CHG
#>   <chr>     <dbl> <dbl> <chr> <dbl> <dbl>
#> 1 1             3   117 NA      113     4
#> 2 3             1   101 NA      111   -10
#> 3 3             2   123 NA      111    12
#> 4 1             0   113 Y       113    NA
#> 5 2             0    95 Y        95    NA
#> 6 3             0   111 Y       111    NA
#> 7 1            -1   113 NA       NA    NA
```
