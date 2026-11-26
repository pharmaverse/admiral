# Get the Key Variables of a Dataset

**\[experimental\]**

Returns the key variables stored for a dataset with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).

## Usage

``` r
get_dataset_keys(dataset)
```

## Arguments

- dataset:

  A data frame

  Permitted values

  :   a data frame

  Default value

  :   none

## Value

The key variables as a list of symbols, in key order, so that they can
be passed to a `by_vars` argument or back to
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md);
an empty list if the dataset has been declared to have no key variables;
or `NULL` if its keys have not been set. The keys are never inferred.

## See also

The admiral_df Class:
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md),
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_class.md),
[`get_dataset_name()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_name.md),
[`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)

## Examples

``` r
library(tibble)

advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "1",      "DIABP",  "BASELINE",    51,
  "1",      "SYSBP",  "BASELINE",   121
)

get_dataset_keys(advs)
#> NULL

advs <- admiral_df(advs, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT))
get_dataset_keys(advs)
#> [[1]]
#> USUBJID
#> 
#> [[2]]
#> PARAMCD
#> 
#> [[3]]
#> AVISIT
#> 
```
