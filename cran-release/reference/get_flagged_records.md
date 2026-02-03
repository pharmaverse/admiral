# Create an Existence Flag

Create a flag variable for the input dataset which indicates if there
exists at least one observation in the input dataset fulfilling a
certain condition.

**Note:** This is a helper function for
`derive_vars_merged_exist_flag()` which inputs this result into
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md).

## Usage

``` r
get_flagged_records(dataset, new_var, condition, filter = NULL)
```

## Arguments

- dataset:

  Input dataset

  Default value

  :   none

- new_var:

  New variable

  The specified variable is added to the input dataset.

  Default value

  :   none

- condition:

  Condition

  The condition is evaluated at the dataset (`dataset`). For all rows
  where it evaluates as `TRUE` the new variable is set to `1` in the new
  column. Otherwise, it is set to `0`.

  Default value

  :   none

- filter:

  Filter for additional data

  Only observations fulfilling the specified condition are taken into
  account for flagging. If the argument is not specified, all
  observations are considered.

  Permitted values

  :   a condition

  Default value

  :   `NULL`

## Value

The output dataset is the input dataset filtered by the `filter`
condition and with the variable specified for `new_var` representing a
flag for the condition.

## See also

Utilities used within Derivation functions:
[`extract_unit()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/extract_unit.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_vars_query.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)


ae <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID,    ~AETERM,     ~AEREL,
  "PILOT01",    "AE", "01-1028", "ERYTHEMA", "POSSIBLE",
  "PILOT01",    "AE", "01-1028", "PRURITUS", "PROBABLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "POSSIBLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "PROBABLE"
)


get_flagged_records(
  dataset = ae,
  new_var = AERELFL,
  condition = AEREL == "PROBABLE"
) %>%
  select(STUDYID, USUBJID, AERELFL)
#> # A tibble: 4 × 3
#>   STUDYID USUBJID AERELFL
#>   <chr>   <chr>     <dbl>
#> 1 PILOT01 01-1028       0
#> 2 PILOT01 01-1028       1
#> 3 PILOT01 06-1049       0
#> 4 PILOT01 06-1049       1

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
get_flagged_records(
  dataset = vs,
  new_var = WTBLHIFL,
  condition = VSSTRESN > 90,
  filter = VSTESTCD == "WEIGHT" & VSBLFL == "Y"
) %>%
  select(STUDYID, USUBJID, WTBLHIFL)
#> # A tibble: 3 × 3
#>   STUDYID USUBJID WTBLHIFL
#>   <chr>   <chr>      <dbl>
#> 1 PILOT01 01-1028        1
#> 2 PILOT01 04-1127        0
#> 3 PILOT01 06-1049        0
```
