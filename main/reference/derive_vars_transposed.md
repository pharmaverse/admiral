# Derive Variables by Transposing and Merging a Second Dataset

Adds variables from a vertical dataset after transposing it into a wide
one.

## Usage

``` r
derive_vars_transposed(
  dataset,
  dataset_merge,
  by_vars,
  id_vars = NULL,
  key_var,
  value_var,
  filter = NULL,
  relationship = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- dataset_merge:

  Dataset to transpose and merge

  The variables specified by the `by_vars`, `id_vars`, `key_var` and
  `value_var` arguments are expected. The variables `by_vars`,
  `id_vars`, `key_var` have to be a unique key.

  Default value

  :   none

- by_vars:

  Grouping variables

  Keys used to merge `dataset_merge` with `dataset`.

  Default value

  :   none

- id_vars:

  ID variables

  Variables (excluding `by_vars` and `key_var`) that uniquely identify
  each observation in `dataset_merge`.

  Default value

  :   `NULL`

- key_var:

  The variable of `dataset_merge` containing the names of the transposed
  variables

  Default value

  :   none

- value_var:

  The variable of `dataset_merge` containing the values of the
  transposed variables

  Default value

  :   none

- filter:

  Expression used to restrict the records of `dataset_merge` prior to
  transposing

  Default value

  :   `NULL`

- relationship:

  Expected merge-relationship between the `by_vars` variable(s) in
  `dataset` and `dataset_merge` (after transposition)

  This argument is passed to the
  [`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  function. See
  <https://dplyr.tidyverse.org/reference/mutate-joins.html#arguments>
  for more details.

  Permitted values

  :   `"one-to-one"`, `"one-to-many"`, `"many-to-one"`,
      `"many-to-many"`, `NULL`

  Default value

  :   `NULL`

## Value

The input dataset with transposed variables from `dataset_merge` added

## Details

1.  The records from the dataset to transpose and merge
    (`dataset_merge`) are restricted to those matching the `filter`
    condition, if provided.

2.  The records from `dataset_merge` are checked to ensure they are
    uniquely identified using `by_vars`, `id_vars` and `key_var`.

3.  `dataset_merge` is transposed (from "tall" to "wide"), with new
    variables added whose names come from `key_var` and values come from
    `value_var`.

4.  The transposed dataset is merged with the input `dataset` using
    `by_vars` as keys. If a `relationship` has been provided, this merge
    must satisfy the relationship, otherwise an error is thrown.

Note that unlike other `derive_vars_*()` functions, the final step may
cause new records to be added to the input dataset. The `relationship`
argument can be specified to ensure this does not happen inadvertently.

## See also

[`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

# Adding ATC classes to CM using FACM
cm <- tribble(
  ~USUBJID,       ~CMGRPID, ~CMREFID,  ~CMDECOD,
  "BP40257-1001", "14",     "1192056", "PARACETAMOL",
  "BP40257-1001", "18",     "2007001", "SOLUMEDROL",
  "BP40257-1002", "19",     "2791596", "SPIRONOLACTONE"
)
facm <- tribble(
  ~USUBJID,       ~FAGRPID, ~FAREFID,  ~FATESTCD,  ~FASTRESC,
  "BP40257-1001", "1",      "1192056", "CMATC1CD", "N",
  "BP40257-1001", "1",      "1192056", "CMATC2CD", "N02",
  "BP40257-1001", "1",      "1192056", "CMATC3CD", "N02B",
  "BP40257-1001", "1",      "1192056", "CMATC4CD", "N02BE",
  "BP40257-1001", "1",      "2007001", "CMATC1CD", "D",
  "BP40257-1001", "1",      "2007001", "CMATC2CD", "D10",
  "BP40257-1001", "1",      "2007001", "CMATC3CD", "D10A",
  "BP40257-1001", "1",      "2007001", "CMATC4CD", "D10AA",
  "BP40257-1001", "2",      "2007001", "CMATC1CD", "D",
  "BP40257-1001", "2",      "2007001", "CMATC2CD", "D07",
  "BP40257-1001", "2",      "2007001", "CMATC3CD", "D07A",
  "BP40257-1001", "2",      "2007001", "CMATC4CD", "D07AA",
  "BP40257-1001", "3",      "2007001", "CMATC1CD", "H",
  "BP40257-1001", "3",      "2007001", "CMATC2CD", "H02",
  "BP40257-1001", "3",      "2007001", "CMATC3CD", "H02A",
  "BP40257-1001", "3",      "2007001", "CMATC4CD", "H02AB",
  "BP40257-1002", "1",      "2791596", "CMATC1CD", "C",
  "BP40257-1002", "1",      "2791596", "CMATC2CD", "C03",
  "BP40257-1002", "1",      "2791596", "CMATC3CD", "C03D",
  "BP40257-1002", "1",      "2791596", "CMATC4CD", "C03DA"
)

cm %>%
  derive_vars_transposed(
    dataset_merge = facm,
    by_vars = exprs(USUBJID, CMREFID = FAREFID),
    id_vars = exprs(FAGRPID),
    key_var = FATESTCD,
    value_var = FASTRESC
  ) %>%
  select(USUBJID, CMDECOD, starts_with("CMATC"))
#> # A tibble: 5 × 6
#>   USUBJID      CMDECOD        CMATC1CD CMATC2CD CMATC3CD CMATC4CD
#>   <chr>        <chr>          <chr>    <chr>    <chr>    <chr>   
#> 1 BP40257-1001 PARACETAMOL    N        N02      N02B     N02BE   
#> 2 BP40257-1001 SOLUMEDROL     D        D10      D10A     D10AA   
#> 3 BP40257-1001 SOLUMEDROL     D        D07      D07A     D07AA   
#> 4 BP40257-1001 SOLUMEDROL     H        H02      H02A     H02AB   
#> 5 BP40257-1002 SPIRONOLACTONE C        C03      C03D     C03DA   

# Note: the `id_vars` argument here is needed to uniquely identify
# rows of dataset_merge and avoid duplicates-related errors.
# Compare the above call to when `id_vars = NULL`:

try(
  cm %>%
    derive_vars_transposed(
      dataset_merge = facm,
      by_vars = exprs(USUBJID, CMREFID = FAREFID),
      id_vars = NULL,
      key_var = FATESTCD,
      value_var = FASTRESC
    )
)
#> Error in signal_duplicate_records(dataset_merge, by_vars = c(by_vars,  : 
#>   Dataset `dataset_merge` contains duplicate records with respect to
#> `USUBJID`, `FAREFID`, and `FATESTCD`
#> Please check data and `by_vars`, `id_vars`, and `key_var` arguments.
#> ℹ Run `admiral::get_duplicates_dataset()` to access the duplicate records
```
