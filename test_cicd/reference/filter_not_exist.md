# Returns records that don't fit into existing by groups in a filtered source dataset

Returns all records in the input dataset that belong to by groups that
are not present in a source dataset, after the source dataset is
optionally filtered. For example, this could be used to return ADSL
records for subjects that didn't take certain concomitant medications
during the course of the study (as per records in ADCM).

## Usage

``` r
filter_not_exist(dataset, dataset_add, by_vars, filter_add = NULL)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- dataset_add:

  Source dataset

  The source dataset, which determines the by groups returned in the
  input dataset, based on the groups that don't exist in this dataset
  after being subset by `filter_add`.

  The variables specified in the `by_vars` and `filter_add` parameters
  are expected in this dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   none

- filter_add:

  Filter for the source dataset

  The filter condition which will be used to subset the source dataset.
  Alternatively, if no filter condition is supplied, no subsetting of
  the source dataset will be performed.

  Default value

  :   `NULL`

## Value

The records in the input dataset which are not contained within any
existing by group in the filtered source dataset.

## Details

Returns the records in `dataset` which don't match any existing by
groups in `dataset_add`, after being filtered according to `filter_add`.
If all by groups that exist in `dataset` don't exist in `dataset_add`,
an empty dataset will be returned.

## See also

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/count_vals.md),
[`filter_exist()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/filter_exist.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/filter_extreme.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/filter_joined.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/min_cond.md)

## Examples

``` r
# Get demographic information about subjects who didn't take vitamin supplements
# during the study

library(tibble)

adsl <- tribble(
  ~USUBJID,      ~AGE, ~SEX,
  "01-701-1015", 63,   "F",
  "01-701-1023", 64,   "M",
  "01-701-1034", 77,   "F",
  "01-701-1118", 52,   "M"
)

adcm <- tribble(
  ~USUBJID,      ~CMTRT,         ~CMSTDTC,
  "01-701-1015", "ASPIRIN",      "2013-05-14",
  "01-701-1023", "MYLANTA",      "2014-01-04",
  "01-701-1023", "CALCIUM",      "2014-02-25",
  "01-701-1034", "VITAMIN C",    "2013-12-12",
  "01-701-1034", "CALCIUM",      "2013-03-27",
  "01-701-1118", "MULTIVITAMIN", "2013-02-21"
)

filter_not_exist(
  dataset = adsl,
  dataset_add = adcm,
  by_vars = exprs(USUBJID),
  filter_add = str_detect(CMTRT, "VITAMIN")
)
#> # A tibble: 2 × 3
#>   USUBJID       AGE SEX  
#>   <chr>       <dbl> <chr>
#> 1 01-701-1015    63 F    
#> 2 01-701-1023    64 M    
```
