# Merge an Existence Flag

Adds a flag variable to the input dataset which indicates if there
exists at least one observation in another dataset fulfilling a certain
condition.

**Note:** This is a wrapper function for the more generic
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md).

## Usage

``` r
derive_var_merged_exist_flag(
  dataset,
  dataset_add,
  by_vars,
  new_var,
  condition,
  true_value = "Y",
  false_value = NA_character_,
  missing_value = NA_character_,
  filter_add = NULL
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

  The variables specified by the `by_vars` argument are expected.

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

- new_var:

  New variable

  The specified variable is added to the input dataset.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   none

- condition:

  Condition

  The condition is evaluated at the additional dataset (`dataset_add`).
  For all by groups where it evaluates as `TRUE` at least once the new
  variable is set to the true value (`true_value`). For all by groups
  where it evaluates as `FALSE` or `NA` for all observations the new
  variable is set to the false value (`false_value`). The new variable
  is set to the missing value (`missing_value`) for by groups not
  present in the additional dataset.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   none

- true_value:

  True value

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"Y"`

- false_value:

  False value

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- missing_value:

  Value used for missing information

  The new variable is set to the specified value for all by groups
  without observations in the additional dataset.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- filter_add:

  Filter for additional data

  Only observations fulfilling the specified condition are taken into
  account for flagging. If the argument is not specified, all
  observations are considered.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variable specified for `new_var` derived
from the additional dataset (`dataset_add`).

## Details

1.  The additional dataset is restricted to the observations matching
    the `filter_add` condition.

2.  The new variable is added to the input dataset and set to the true
    value (`true_value`) if for the by group at least one observation
    exists in the (restricted) additional dataset where the condition
    evaluates to `TRUE`. It is set to the false value (`false_value`) if
    for the by group at least one observation exists and for all
    observations the condition evaluates to `FALSE` or `NA`. Otherwise,
    it is set to the missing value (`missing_value`).

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_ef_msrc.md),
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
library(dplyr, warn.conflicts = FALSE)

dm <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1028",   71, "YEARS",
  "PILOT01",    "DM", "04-1127",   84, "YEARS",
  "PILOT01",    "DM", "06-1049",   60, "YEARS"
)

ae <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID,    ~AETERM,     ~AEREL,
  "PILOT01",    "AE", "01-1028", "ERYTHEMA", "POSSIBLE",
  "PILOT01",    "AE", "01-1028", "PRURITUS", "PROBABLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "POSSIBLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "PROBABLE"
)


derive_var_merged_exist_flag(
  dm,
  dataset_add = ae,
  by_vars = exprs(STUDYID, USUBJID),
  new_var = AERELFL,
  condition = AEREL == "PROBABLE"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, AERELFL)
#> # A tibble: 3 × 5
#>   STUDYID USUBJID   AGE AGEU  AERELFL
#>   <chr>   <chr>   <dbl> <chr> <chr>  
#> 1 PILOT01 01-1028    71 YEARS Y      
#> 2 PILOT01 04-1127    84 YEARS NA     
#> 3 PILOT01 06-1049    60 YEARS Y      

vs <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSBLFL,
  "PILOT01",    "VS", "01-1028", "SCREENING",  "HEIGHT",     177.8,      NA,
  "PILOT01",    "VS", "01-1028", "SCREENING",  "WEIGHT",     98.88,      NA,
  "PILOT01",    "VS", "01-1028",  "BASELINE",  "WEIGHT",     99.34,     "Y",
  "PILOT01",    "VS", "01-1028",    "WEEK 4",  "WEIGHT",     98.88,      NA,
  "PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,      NA,
  "PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,      NA,
  "PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,     "Y",
  "PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,      NA,
  "PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,      NA,
  "PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,      NA,
  "PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     "Y",
  "PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,      NA
)
derive_var_merged_exist_flag(
  dm,
  dataset_add = vs,
  by_vars = exprs(STUDYID, USUBJID),
  filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
  new_var = WTBLHIFL,
  condition = VSSTRESN > 90,
  false_value = "N",
  missing_value = "M"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, WTBLHIFL)
#> # A tibble: 3 × 5
#>   STUDYID USUBJID   AGE AGEU  WTBLHIFL
#>   <chr>   <chr>   <dbl> <chr> <chr>   
#> 1 PILOT01 01-1028    71 YEARS Y       
#> 2 PILOT01 04-1127    84 YEARS N       
#> 3 PILOT01 06-1049    60 YEARS N       
```
