# Adds a Variable Numbering the Observations Within Each By Group

Adds a variable numbering the observations within each by group

## Usage

``` r
derive_var_obs_number(
  dataset,
  by_vars = NULL,
  order = NULL,
  new_var = ASEQ,
  check_type = "none"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   `NULL`

- order:

  Sort order

  Within each by group the observations are ordered by the specified
  order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/test_cicd/articles/generic.md).

  Permitted values

  :   list of variables or functions of variables

  Default value

  :   `NULL`

- new_var:

  Name of variable to create

  The new variable is set to the observation number for each by group.
  The numbering starts with 1.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `ASEQ`

- check_type:

  Check uniqueness?

  If `"message"`, `"warning"` or `"error"` is specified, the specified
  message is issued if the observations of the input dataset are not
  unique with respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"none"`

## Value

A dataset containing all observations and variables of the input dataset
and additionally the variable specified by the `new_var` parameter.

## Details

For each group (with respect to the variables specified for the
`by_vars` parameter) the first or last observation (with respect to the
order specified for the `order` parameter and the mode specified for the
`mode` parameter) is included in the output dataset.

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_exist_flag.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_transposed.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)
vs <- tribble(
  ~STUDYID,  ~DOMAIN,      ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
  "PILOT01",    "VS", "01-703-1182",   "DIABP",         3,       815,
  "PILOT01",    "VS", "01-703-1182",   "DIABP",         3,       816,
  "PILOT01",    "VS", "01-703-1182",   "DIABP",         4,       815,
  "PILOT01",    "VS", "01-703-1182",   "DIABP",         4,       816,
  "PILOT01",    "VS", "01-703-1182",   "PULSE",         3,       815,
  "PILOT01",    "VS", "01-703-1182",   "PULSE",         3,       816,
  "PILOT01",    "VS", "01-703-1182",   "PULSE",         4,       815,
  "PILOT01",    "VS", "01-703-1182",   "PULSE",         4,       816,
  "PILOT01",    "VS", "01-703-1182",   "SYSBP",         3,       815,
  "PILOT01",    "VS", "01-703-1182",   "SYSBP",         3,       816,
  "PILOT01",    "VS", "01-703-1182",   "SYSBP",         4,       815,
  "PILOT01",    "VS", "01-703-1182",   "SYSBP",         4,       816,
  "PILOT01",    "VS", "01-716-1229",   "DIABP",         3,       815,
  "PILOT01",    "VS", "01-716-1229",   "DIABP",         3,       816,
  "PILOT01",    "VS", "01-716-1229",   "DIABP",         4,       815,
  "PILOT01",    "VS", "01-716-1229",   "DIABP",         4,       816,
  "PILOT01",    "VS", "01-716-1229",   "PULSE",         3,       815,
  "PILOT01",    "VS", "01-716-1229",   "PULSE",         3,       816,
  "PILOT01",    "VS", "01-716-1229",   "PULSE",         4,       815,
  "PILOT01",    "VS", "01-716-1229",   "PULSE",         4,       816,
  "PILOT01",    "VS", "01-716-1229",   "SYSBP",         3,       815,
  "PILOT01",    "VS", "01-716-1229",   "SYSBP",         3,       816,
  "PILOT01",    "VS", "01-716-1229",   "SYSBP",         4,       815,
  "PILOT01",    "VS", "01-716-1229",   "SYSBP",         4,       816
)
vs %>%
  derive_var_obs_number(
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VISITNUM, desc(VSTPTNUM))
  )
#> # A tibble: 24 × 7
#>    STUDYID DOMAIN USUBJID     VSTESTCD VISITNUM VSTPTNUM  ASEQ
#>    <chr>   <chr>  <chr>       <chr>       <dbl>    <dbl> <int>
#>  1 PILOT01 VS     01-703-1182 DIABP           3      816     1
#>  2 PILOT01 VS     01-703-1182 DIABP           3      815     2
#>  3 PILOT01 VS     01-703-1182 DIABP           4      816     3
#>  4 PILOT01 VS     01-703-1182 DIABP           4      815     4
#>  5 PILOT01 VS     01-703-1182 PULSE           3      816     1
#>  6 PILOT01 VS     01-703-1182 PULSE           3      815     2
#>  7 PILOT01 VS     01-703-1182 PULSE           4      816     3
#>  8 PILOT01 VS     01-703-1182 PULSE           4      815     4
#>  9 PILOT01 VS     01-703-1182 SYSBP           3      816     1
#> 10 PILOT01 VS     01-703-1182 SYSBP           3      815     2
#> # ℹ 14 more rows
```
