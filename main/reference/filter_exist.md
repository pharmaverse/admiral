# Returns records that fit into existing by groups in a filtered source dataset

Returns all records in the input dataset that belong to by groups that
are present in a source dataset, after the source dataset is optionally
filtered. For example, this could be used to return ADSL records for
subjects that experienced a certain adverse event during the course of
the study (as per records in ADAE).

## Usage

``` r
filter_exist(dataset, dataset_add, by_vars, filter_add = NULL)
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
  input dataset, based on the groups that exist in this dataset after
  being subset by `filter_add`.

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

The records in the input dataset which are contained within an existing
by group in the filtered source dataset.

## Details

Returns the records in `dataset` which match an existing by group in
`dataset_add`, after being filtered according to `filter_add`. If there
are no by groups that exist in both datasets, an empty dataset will be
returned.

## See also

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/main/reference/filter_extreme.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_not_exist.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/main/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)

## Examples

``` r
# Get demographic information about subjects who have suffered from moderate or
# severe fatigue

library(tibble)

adsl <- tribble(
  ~USUBJID,      ~AGE, ~SEX,
  "01-701-1015", 63,   "F",
  "01-701-1034", 77,   "F",
  "01-701-1115", 84,   "M",
  "01-701-1146", 75,   "F",
  "01-701-1444", 63,   "M"
)

adae <- tribble(
  ~USUBJID,      ~AEDECOD,                    ~AESEV,     ~AESTDTC,
  "01-701-1015", "DIARRHOEA",                 "MODERATE", "2014-01-09",
  "01-701-1034", "FATIGUE",                   "SEVERE",   "2014-11-02",
  "01-701-1034", "APPLICATION SITE PRURITUS", "MODERATE", "2014-08-27",
  "01-701-1115", "FATIGUE",                   "MILD",     "2013-01-14",
  "01-701-1146", "FATIGUE",                   "MODERATE", "2013-06-03"
)

filter_exist(
  dataset = adsl,
  dataset_add = adae,
  by_vars = exprs(USUBJID),
  filter_add = AEDECOD == "FATIGUE" & AESEV %in% c("MODERATE", "SEVERE")
)
#> # A tibble: 2 × 3
#>   USUBJID       AGE SEX  
#>   <chr>       <dbl> <chr>
#> 1 01-701-1034    77 F    
#> 2 01-701-1146    75 F    
```
