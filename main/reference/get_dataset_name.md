# Get the Name of a Dataset

**\[experimental\]**

Returns the dataset name stored for a dataset with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md).

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

Utilities used within Derivation functions:
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md),
[`extract_unit()`](https:/pharmaverse.github.io/admiral/main/reference/extract_unit.md),
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/get_dataset_class.md),
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/main/reference/get_dataset_keys.md),
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/main/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/get_vars_query.md)

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
