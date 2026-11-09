# Get the ADaM Dataset Class

**\[experimental\]**

Returns the ADaM class declared for a dataset with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md)
if there is one, and otherwise the class inferred from its variables and
records, see
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md).

## Usage

``` r
get_dataset_class(dataset)
```

## Arguments

- dataset:

  A data frame

  Permitted values

  :   a data frame

  Default value

  :   none

## Value

A character scalar: the declared class, or the inferred one (which may
be `"other"` or `NA`, see
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md)).

## See also

Utilities used within Derivation functions:
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md),
[`extract_unit()`](https:/pharmaverse.github.io/admiral/main/reference/extract_unit.md),
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/main/reference/get_dataset_keys.md),
[`get_dataset_name()`](https:/pharmaverse.github.io/admiral/main/reference/get_dataset_name.md),
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/main/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/get_vars_query.md)

## Examples

``` r
library(tibble)

advs <- tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL,
  "P",      "1",      "DIABP",     51,
  "P",      "1",      "SYSBP",    121
)

# no class has been declared, so it is inferred
get_dataset_class(advs)
#> [1] "BDS"

# a declared class takes precedence
advs <- admiral_df(advs, dataset_class = "TTE")
get_dataset_class(advs)
#> [1] "TTE"
```
