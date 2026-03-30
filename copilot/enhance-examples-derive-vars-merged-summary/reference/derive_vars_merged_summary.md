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
      [`exprs()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/reexport-exprs.md),
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/reexport-exprs.md),
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/reexport-exprs.md),
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

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_summary_records.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/get_summary_records.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_merged_lookup.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/copilot/enhance-examples-derive-vars-merged-summary/reference/derive_vars_transposed.md)

## Examples

### Data setup

The following examples use the BDS dataset below as a basis.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)

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

### Summarize one or more variables using summary functions (`new_vars`)

The `new_vars` argument specifies a named list of expressions where the
right-hand side uses summary functions (e.g.
[`mean()`](https://rdrr.io/r/base/mean.html),
[`sum()`](https://rdrr.io/r/base/sum.html),
[`max()`](https://rdrr.io/r/base/Extremes.html)) to aggregate values
from `dataset_add` within each by group. Multiple summary variables can
be added in a single call.

In the example below, the mean and sum of `AVAL` within each subject and
visit are derived and merged back onto the input dataset:

    derive_vars_merged_summary(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(USUBJID, AVISIT),
      new_vars = exprs(
        MEANVIS = mean(AVAL, na.rm = TRUE),
        SUMVIS = sum(AVAL, na.rm = TRUE)
      )
    )
    #> # A tibble: 9 × 6
    #>   USUBJID AVISIT  ASEQ  AVAL MEANVIS SUMVIS
    #>   <chr>   <chr>  <dbl> <dbl>   <dbl>  <dbl>
    #> 1 1       WEEK 1     1    10      10     10
    #> 2 1       WEEK 1     2    NA      10     10
    #> 3 1       WEEK 2     3    NA     NaN      0
    #> 4 1       WEEK 3     4    42      42     42
    #> 5 1       WEEK 4     5    12      13     39
    #> 6 1       WEEK 4     6    12      13     39
    #> 7 1       WEEK 4     7    15      13     39
    #> 8 2       WEEK 1     1    21      21     21
    #> 9 2       WEEK 4     2    22      22     22

In the example above, subject `"1"` at `"WEEK 2"` has only missing
`AVAL` values, so `MEANVIS` is `NaN` (the result of
`mean(NA, na.rm = TRUE)`) and `SUMVIS` is `0`.

### Restricting source records (`filter_add`)

The `filter_add` argument restricts the records in `dataset_add` that
are used for the summarization. Only records satisfying the filter
condition contribute to the summary values. This can be useful, for
example, to compute a summary statistic based only on records before or
after a certain time point.

In the following example, the mean of `AVAL` is computed only for
records with a positive study day (`ADY > 0`), and the result is merged
onto the `ADSL`-like dataset. Subject `"2"` has no records with
`ADY > 0`, so `MEANPBL` is `NA` for that subject.

    adsl <- tribble(
      ~USUBJID,
      "1",
      "2",
      "3"
    )

    adbds2 <- tribble(
      ~USUBJID, ~ADY, ~AVAL,
      "1",        -3,    10,
      "1",         2,    12,
      "1",         8,    15,
      "3",         4,    42
    )

    derive_vars_merged_summary(
      adsl,
      dataset_add = adbds2,
      by_vars = exprs(USUBJID),
      new_vars = exprs(MEANPBL = mean(AVAL, na.rm = TRUE)),
      filter_add = ADY > 0
    )
    #> # A tibble: 3 × 2
    #>   USUBJID MEANPBL
    #>   <chr>     <dbl>
    #> 1 1          13.5
    #> 2 2          NA
    #> 3 3          42  

### Handling non-matching observations (`missing_values`)

By default, records in `dataset` with no matching by group in
`dataset_add` receive `NA` for the new variables. The `missing_values`
argument allows you to specify a different value for these non-matching
records.

A common use-case is **population median imputation**: subjects without
a baseline measurement are assigned the median baseline value observed
across the rest of the population. In the example below, `adlb_bl`
contains one baseline record per subject for subjects 1–10, except
subjects `"5"` and `"8"` who have no baseline record. Without
`missing_values` those two subjects would receive `NA`; supplying the
pre-computed population median imputes that value instead:

    adsl2 <- tribble(
      ~USUBJID,
      "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"
    )

    adlb_bl <- tribble(
      ~USUBJID, ~ABLFL, ~AVAL,
      "1",      "Y",       10,
      "2",      "Y",       15,
      "3",      "Y",       20,
      "4",      "Y",       28,
      "6",      "Y",       35,
      "7",      "Y",       42,
      "9",      "Y",       50,
      "10",     "Y",       60
    )

    pop_median <- median(adlb_bl$AVAL, na.rm = TRUE)

    derive_vars_merged_summary(
      adsl2,
      dataset_add = adlb_bl,
      by_vars = exprs(USUBJID),
      new_vars = exprs(BASE = mean(AVAL, na.rm = TRUE)),
      filter_add = ABLFL == "Y",
      missing_values = exprs(BASE = !!pop_median)
    )
    #> # A tibble: 10 × 2
    #>    USUBJID  BASE
    #>    <chr>   <dbl>
    #>  1 1        10
    #>  2 2        15
    #>  3 3        20
    #>  4 4        28
    #>  5 5        31.5
    #>  6 6        35
    #>  7 7        42
    #>  8 8        31.5
    #>  9 9        50
    #> 10 10       60  

### Renaming by variables (`by_vars`)

The `by_vars` argument supports renaming, using the syntax
`exprs(<left_name> = <right_name>)`, where `<left_name>` is the variable
name in `dataset` and `<right_name>` is the corresponding variable in
`dataset_add`. This is useful when the grouping variable has different
names in the two datasets.

In the example below the input dataset uses `AVISIT` while the
additional dataset uses `VISIT` for the same concept. The `by_vars`
argument maps them together so the merge can proceed correctly:

    adbds_renamed <- adbds %>% rename(VISIT = AVISIT)

    derive_vars_merged_summary(
      adbds,
      dataset_add = adbds_renamed,
      by_vars = exprs(USUBJID, AVISIT = VISIT),
      new_vars = exprs(MEANVIS = mean(AVAL, na.rm = TRUE))
    )
    #> # A tibble: 9 × 5
    #>   USUBJID AVISIT  ASEQ  AVAL MEANVIS
    #>   <chr>   <chr>  <dbl> <dbl>   <dbl>
    #> 1 1       WEEK 1     1    10      10
    #> 2 1       WEEK 1     2    NA      10
    #> 3 1       WEEK 2     3    NA     NaN
    #> 4 1       WEEK 3     4    42      42
    #> 5 1       WEEK 4     5    12      13
    #> 6 1       WEEK 4     6    12      13
    #> 7 1       WEEK 4     7    15      13
    #> 8 2       WEEK 1     1    21      21
    #> 9 2       WEEK 4     2    22      22

### String aggregation

Summary expressions are not restricted to numeric aggregations. Any
expression that reduces a group to a single value is permitted. For
example, `paste(..., collapse = ", ")` can be used to concatenate
character values within a by group into a single string.

In the example below, the lesion identifiers observed at baseline for
each subject are collected into a single comma-separated string and
merged onto the `ADSL` dataset:

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
    #> 3 3       <NA>                          
