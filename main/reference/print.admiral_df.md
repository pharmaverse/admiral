# Print `admiral_df` Objects

**\[experimental\]**

Prints the dataset name, ADaM class, and key variables stored by
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md)
in a header line, followed by the dataset as it would print without the
`admiral_df` class. The header is only printed if any of them have been
stored, so that the output of a derivation which merely carries the
class prints as before.

If no class has been stored, the class inferred from the dataset (see
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md))
is shown instead, marked as inferred, provided it is one of the ADaM
dataset classes.

## Usage

``` r
# S3 method for class 'admiral_df'
print(x, ...)
```

## Arguments

- x:

  An `admiral_df` object

  Default value

  :   none

- ...:

  Passed on to the next print method

  Default value

  :   none

## Value

`x`, invisibly

## See also

Utilities for printing:
[`print.adam_templates()`](https:/pharmaverse.github.io/admiral/main/reference/print.adam_templates.md),
[`print.duplicates()`](https:/pharmaverse.github.io/admiral/main/reference/print.duplicates.md),
[`print.source()`](https:/pharmaverse.github.io/admiral/main/reference/print.source.md),
[`print_named_list()`](https:/pharmaverse.github.io/admiral/main/reference/print_named_list.md)

## Examples

``` r
library(tibble)

advs <- tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "P",      "1",      "DIABP",  "BASELINE",    51,
  "P",      "1",      "SYSBP",  "BASELINE",   121,
  "P",      "2",      "DIABP",  "BASELINE",    79
)

admiral_df(
  advs,
  dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
  dataset_name = "ADVS",
  dataset_class = "BDS"
)
#> # An admiral dataset (name: ADVS, class: BDS, keys: USUBJID, PARAMCD, AVISIT)
#> # A tibble: 3 × 5
#>   STUDYID USUBJID PARAMCD AVISIT    AVAL
#>   <chr>   <chr>   <chr>   <chr>    <dbl>
#> 1 P       1       DIABP   BASELINE    51
#> 2 P       1       SYSBP   BASELINE   121
#> 3 P       2       DIABP   BASELINE    79
```
