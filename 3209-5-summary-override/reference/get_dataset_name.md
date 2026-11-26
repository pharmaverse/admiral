# Get the Name of a Dataset

**\[experimental\]**

Returns the dataset name stored for a dataset with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).

## Usage

``` r
get_dataset_name(dataset)
```

## Arguments

- dataset:

  A data frame

  Permitted values

  :   a data frame

  Default value

  :   none

## Value

A character scalar, or `NULL` if no name has been set.

## See also

The admiral_df Class:
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md),
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_class.md),
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_keys.md),
[`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)

## Examples

``` r
library(tibble)

advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL,
  "1",      "DIABP",     51
)

get_dataset_name(advs)
#> NULL

advs <- admiral_df(advs, dataset_name = "ADVS")
get_dataset_name(advs)
#> [1] "ADVS"
```
