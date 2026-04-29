# Merge an Existence Flag From Multiple Sources

Adds a flag variable to the input dataset which indicates if there
exists at least one observation in one of the source datasets fulfilling
a certain condition. For example, if a dose adjustment flag should be
added to `ADEX` but the dose adjustment information is collected in
different datasets, e.g., `EX`, `EC`, and `FA`.

## Usage

``` r
derive_var_merged_ef_msrc(
  dataset,
  by_vars,
  flag_events,
  source_datasets,
  new_var,
  true_value = "Y",
  false_value = NA_character_,
  missing_value = NA_character_
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

- by_vars:

  Grouping variables

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- flag_events:

  Flag events

  A list of
  [`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md)
  objects is expected. For each event the condition (`condition` field)
  is evaluated in the source dataset referenced by the `dataset_name`
  field. If it evaluates to `TRUE` at least once, the new variable is
  set to `true_value`.

  Permitted values

  :   a list of
      [`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md)
      objects

  Default value

  :   none

- source_datasets:

  Source datasets

  A named list of datasets is expected. The `dataset_name` field of
  [`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md)
  refers to the dataset provided in the list.

  Permitted values

  :   named list of datasets, e.g., `list(adsl = adsl, ae = ae)`

  Default value

  :   none

- new_var:

  New variable

  The specified variable is added to the input dataset.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   none

- true_value:

  True value

  The new variable (`new_var`) is set to the specified value for all by
  groups for which at least one of the source object (`sources`) has the
  condition evaluate to `TRUE`.

  The values of `true_value`, `false_value`, and `missing_value` must be
  of the same type.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"Y"`

- false_value:

  False value

  The new variable (`new_var`) is set to the specified value for all by
  groups which occur in at least one source (`sources`) but the
  condition never evaluates to `TRUE`.

  The values of `true_value`, `false_value`, and `missing_value` must be
  of the same type.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- missing_value:

  Values used for missing information

  The new variable is set to the specified value for all by groups
  without observations in any of the sources (`sources`).

  The values of `true_value`, `false_value`, and `missing_value` must be
  of the same type.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variable specified for `new_var`.

## Details

1.  For each
    [`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md)
    object specified for `flag_events`: The condition (`condition`) is
    evaluated in the dataset referenced by `dataset_name`. If the
    `by_vars` field is specified the dataset is grouped by the specified
    variables for evaluating the condition. If named elements are used
    in `by_vars` like `by_vars = exprs(USUBJID, EXLNKID = ECLNKID)`, the
    variables are renamed after the evaluation. If the `by_vars` element
    is not specified, the observations are grouped by the variables
    specified for the `by_vars` argument.

2.  The new variable (`new_var`) is added to the input dataset and set
    to the true value (`true_value`) if for the by group at least one
    condition evaluates to `TRUE` in one of the sources. It is set to
    the false value (`false_value`) if for the by group at least one
    observation exists and for all observations the condition evaluates
    to `FALSE` or `NA`. Otherwise, it is set to the missing value
    (`missing_value`).

## See also

[`flag_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/flag_event.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_transposed.md)

## Examples

``` r
library(dplyr)

# Derive a flag indicating anti-cancer treatment based on CM and PR
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3",
  "4"
)

cm <- tribble(
  ~USUBJID, ~CMCAT,        ~CMSEQ,
  "1",      "ANTI-CANCER",      1,
  "1",      "GENERAL",          2,
  "2",      "GENERAL",          1,
  "3",      "ANTI-CANCER",      1
)

# Assuming all records in PR indicate cancer treatment
pr <- tibble::tribble(
  ~USUBJID, ~PRSEQ,
  "2",      1,
  "3",      1
)

derive_var_merged_ef_msrc(
  adsl,
  by_vars = exprs(USUBJID),
  flag_events = list(
    flag_event(
      dataset_name = "cm",
      condition = CMCAT == "ANTI-CANCER"
    ),
    flag_event(
      dataset_name = "pr"
    )
  ),
  source_datasets = list(cm = cm, pr = pr),
  new_var = CANCTRFL
)
#> # A tibble: 4 × 2
#>   USUBJID CANCTRFL
#>   <chr>   <chr>   
#> 1 1       Y       
#> 2 2       Y       
#> 3 3       Y       
#> 4 4       NA      

# Using different by variables depending on the source
# Add a dose adjustment flag to ADEX based on ADEX, EC, and FA
adex <- tribble(
  ~USUBJID, ~EXLNKID, ~EXADJ,
  "1",      "1",      "AE",
  "1",      "2",      NA_character_,
  "1",      "3",      NA_character_,
  "2",      "1",      NA_character_,
  "3",      "1",      NA_character_
)

ec <- tribble(
  ~USUBJID, ~ECLNKID, ~ECADJ,
  "1",      "3",      "AE",
  "3",      "1",      NA_character_
)

fa <- tribble(
  ~USUBJID, ~FALNKID, ~FATESTCD, ~FAOBJ,            ~FASTRESC,
  "3",      "1",      "OCCUR",   "DOSE ADJUSTMENT", "Y"
)

derive_var_merged_ef_msrc(
  adex,
  by_vars = exprs(USUBJID, EXLNKID),
  flag_events = list(
    flag_event(
      dataset_name = "ex",
      condition = !is.na(EXADJ)
    ),
    flag_event(
      dataset_name = "ec",
      condition = !is.na(ECADJ),
      by_vars = exprs(USUBJID, EXLNKID = ECLNKID)
    ),
    flag_event(
      dataset_name = "fa",
      condition = FATESTCD == "OCCUR" & FAOBJ == "DOSE ADJUSTMENT" & FASTRESC == "Y",
      by_vars = exprs(USUBJID, EXLNKID = FALNKID)
    )
  ),
  source_datasets = list(ex = adex, ec = ec, fa = fa),
  new_var = DOSADJFL
)
#> # A tibble: 5 × 4
#>   USUBJID EXLNKID EXADJ DOSADJFL
#>   <chr>   <chr>   <chr> <chr>   
#> 1 1       1       AE    Y       
#> 2 1       2       NA    NA      
#> 3 1       3       NA    Y       
#> 4 2       1       NA    NA      
#> 5 3       1       NA    Y       
```
