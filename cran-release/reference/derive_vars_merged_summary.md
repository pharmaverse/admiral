# Merge Summary Variables

Merge a summary variable from a dataset to the input dataset.

## Usage

``` r
derive_vars_merged_summary(
  dataset,
  dataset_add,
  by_vars,
  new_vars = NULL,
  filter_add = NULL,
  missing_values = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_add:

  Additional dataset

  The variables specified by the `by_vars` and the variables used on the
  left hand sides of the `new_vars` arguments are expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  The expressions on the left hand sides of `new_vars` are evaluated by
  the specified *variables*. Then the resulting values are merged to the
  input dataset (`dataset`) by the specified *variables*.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- new_vars:

  New variables to add

  The specified variables are added to the input dataset.

  A named list of expressions is expected:

  - LHS refer to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, an expression or NA. If summary
    functions are used, the values are summarized by the variables
    specified for `by_vars`.

  For example:

        new_vars = exprs(
          DOSESUM = sum(AVAL),
          DOSEMEAN = mean(AVAL)
        )

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations fulfilling the specified condition are taken into
  account for summarizing. If the argument is not specified, all
  observations are considered.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- missing_values:

  Values for non-matching observations

  For observations of the input dataset (`dataset`) which do not have a
  matching observation in the additional dataset (`dataset_add`) the
  values of the specified variables are set to the specified value. Only
  variables specified for `new_vars` can be specified for
  `missing_values`.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   `NULL`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variables specified for `new_vars`.

## Details

1.  The records from the additional dataset (`dataset_add`) are
    restricted to those matching the `filter_add` condition.

2.  The new variables (`new_vars`) are created for each by group
    (`by_vars`) in the additional dataset (`dataset_add`) by calling
    [`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html).
    I.e., all observations of a by group are summarized to a single
    observation.

3.  The new variables are merged to the input dataset. For observations
    without a matching observation in the additional dataset the new
    variables are set to `NA`. Observations in the additional dataset
    which have no matching observation in the input dataset are ignored.

## See also

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_summary_records.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_summary_records.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_lookup.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_transposed.md)

## Examples

``` r
library(tibble)

# Add a variable for the mean of AVAL within each visit
adbds <- tribble(
  ~USUBJID,  ~AVISIT,  ~ASEQ, ~AVAL,
  "1",      "WEEK 1",      1,    10,
  "1",      "WEEK 1",      2,    NA,
  "1",      "WEEK 2",      3,    NA,
  "1",      "WEEK 3",      4,    42,
  "1",      "WEEK 4",      5,    12,
  "1",      "WEEK 4",      6,    12,
  "1",      "WEEK 4",      7,    15,
  "2",      "WEEK 1",      1,    21,
  "2",      "WEEK 4",      2,    22
)

derive_vars_merged_summary(
  adbds,
  dataset_add = adbds,
  by_vars = exprs(USUBJID, AVISIT),
  new_vars = exprs(
    MEANVIS = mean(AVAL, na.rm = TRUE),
    MAXVIS = max(AVAL, na.rm = TRUE)
  )
)
#> Warning: There was 1 warning in `reframe()`.
#> ℹ In argument: `MAXVIS = max(AVAL, na.rm = TRUE)`.
#> ℹ In group 2: `USUBJID = "1"` `AVISIT = "WEEK 2"`.
#> Caused by warning in `max()`:
#> ! no non-missing arguments to max; returning -Inf
#> # A tibble: 9 × 6
#>   USUBJID AVISIT  ASEQ  AVAL MEANVIS MAXVIS
#>   <chr>   <chr>  <dbl> <dbl>   <dbl>  <dbl>
#> 1 1       WEEK 1     1    10      10     10
#> 2 1       WEEK 1     2    NA      10     10
#> 3 1       WEEK 2     3    NA     NaN   -Inf
#> 4 1       WEEK 3     4    42      42     42
#> 5 1       WEEK 4     5    12      13     15
#> 6 1       WEEK 4     6    12      13     15
#> 7 1       WEEK 4     7    15      13     15
#> 8 2       WEEK 1     1    21      21     21
#> 9 2       WEEK 4     2    22      22     22

# Add a variable listing the lesion ids at baseline
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
)

adtr <- tribble(
  ~USUBJID,     ~AVISIT, ~LESIONID,
  "1",       "BASELINE",  "INV-T1",
  "1",       "BASELINE",  "INV-T2",
  "1",       "BASELINE",  "INV-T3",
  "1",       "BASELINE",  "INV-T4",
  "1",         "WEEK 1",  "INV-T1",
  "1",         "WEEK 1",  "INV-T2",
  "1",         "WEEK 1",  "INV-T4",
  "2",       "BASELINE",  "INV-T1",
  "2",       "BASELINE",  "INV-T2",
  "2",       "BASELINE",  "INV-T3",
  "2",         "WEEK 1",  "INV-T1",
  "2",         "WEEK 1",  "INV-N1"
)

derive_vars_merged_summary(
  adsl,
  dataset_add = adtr,
  by_vars = exprs(USUBJID),
  filter_add = AVISIT == "BASELINE",
  new_vars = exprs(LESIONSBL = paste(LESIONID, collapse = ", "))
)
#> # A tibble: 3 × 2
#>   USUBJID LESIONSBL                     
#>   <chr>   <chr>                         
#> 1 1       INV-T1, INV-T2, INV-T3, INV-T4
#> 2 2       INV-T1, INV-T2, INV-T3        
#> 3 3       NA                            
```
