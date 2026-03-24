# Execute a Derivation with Different Arguments for Subsets of the Input Dataset

The input dataset is split into slices (subsets) and for each slice the
derivation is called separately. Some or all arguments of the derivation
may vary depending on the slice.

## Usage

``` r
slice_derivation(dataset, derivation, ..., args = NULL)
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
  object passed to the `args` argument or be provided in *every*
  [`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md).

  Default value

  :   none

- ...:

  A
  [`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
  object is expected

  Each slice defines a subset of the input dataset and some of the
  parameters for the derivation. The derivation is called on the subset
  with the parameters specified by the `args` parameter and the `args`
  field of the
  [`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
  object. If a parameter is specified for both, the value in
  [`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
  overwrites the one in `args`.

  Default value

  :   none

- args:

  Arguments of the derivation

  A
  [`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
  object is expected.

  Default value

  :   `NULL`

## Value

The input dataset with the variables derived by the derivation added

## Details

For each slice the derivation is called on the subset defined by the
`filter` field of the
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
object and with the parameters specified by the `args` parameter and the
`args` field of the
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
object. If a parameter is specified for both, the value in
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md)
overwrites the one in `args`.

- Observations that match with more than one slice are only considered
  for the first matching slice.

- The derivation is called for slices with no observations.

- Observations with no match to any of the slices are included in the
  output dataset but the derivation is not called for them.

It is also possible to pass functions from outside the `{admiral}`
package to `slice_derivation()`, e.g. an extension package function, or
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
The only requirement for a function being passed to `derivation` is that
it must take a dataset as its first argument and return a dataset.

## See also

[`params()`](https:/pharmaverse.github.io/admiral/main/reference/params.md)
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)
[`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md)

Higher Order Functions:
[`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md),
[`derivation_slice()`](https:/pharmaverse.github.io/admiral/main/reference/derivation_slice.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)

## Examples

``` r
library(tibble)
library(stringr)
advs <- tribble(
  ~USUBJID, ~VSDTC,       ~VSTPT,
  "1",      "2020-04-16", NA_character_,
  "1",      "2020-04-16", "BEFORE TREATMENT"
)

# For the second slice filter is set to TRUE. Thus derive_vars_dtm is called
# with time_imputation = "last" for all observations which do not match for the
# first slice.
slice_derivation(
  advs,
  derivation = derive_vars_dtm,
  args = params(
    dtc = VSDTC,
    new_vars_prefix = "A"
  ),
  derivation_slice(
    filter = str_detect(VSTPT, "PRE|BEFORE"),
    args = params(time_imputation = "first")
  ),
  derivation_slice(
    filter = TRUE,
    args = params(time_imputation = "last")
  )
)
#> # A tibble: 2 × 5
#>   USUBJID VSDTC      VSTPT            ADTM                ATMF 
#>   <chr>   <chr>      <chr>            <dttm>              <chr>
#> 1 1       2020-04-16 BEFORE TREATMENT 2020-04-16 00:00:00 H    
#> 2 1       2020-04-16 NA               2020-04-16 23:59:59 H    
```
