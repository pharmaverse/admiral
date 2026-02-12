# Derive Baseline Variables

Derive baseline variables, e.g. `BASE` or `BNRIND`, in a BDS dataset.

**Note:** This is a wrapper function for the more generic
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md).

## Usage

``` r
derive_var_base(
  dataset,
  by_vars,
  source_var = AVAL,
  new_var = BASE,
  filter = ABLFL == "Y"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `source_var` arguments
  are expected to be in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Grouping variables uniquely identifying a set of records for which to
  calculate `new_var`.

  Default value

  :   none

- source_var:

  The column from which to extract the baseline value, e.g. `AVAL`

  Default value

  :   `AVAL`

- new_var:

  The name of the newly created baseline column, e.g. `BASE`

  Default value

  :   `BASE`

- filter:

  The condition used to filter `dataset` for baseline records.

  By default `ABLFL == "Y"`

  Default value

  :   `ABLFL == "Y"`

## Value

A new `data.frame` containing all records and variables of the input
dataset plus the `new_var` variable

## Details

For each `by_vars` group, the baseline record is identified by the
condition specified in `filter` which defaults to `ABLFL == "Y"`.
Subsequently, every value of the `new_var` variable for the `by_vars`
group is set to the value of the `source_var` variable of the baseline
record. In case there are multiple baseline records within `by_vars` an
error is issued.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr_dir.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)

dataset <- tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC,   ~AVISIT,    ~ABLFL, ~ANRIND,
  "TEST01", "PAT01", "PARAM01", 10.12, NA,       "Baseline", "Y",    "NORMAL",
  "TEST01", "PAT01", "PARAM01", 9.700, NA,       "Day 7",    NA,     "LOW",
  "TEST01", "PAT01", "PARAM01", 15.01, NA,       "Day 14",   NA,     "HIGH",
  "TEST01", "PAT01", "PARAM02", 8.350, NA,       "Baseline", "Y",    "LOW",
  "TEST01", "PAT01", "PARAM02",    NA, NA,       "Day 7",    NA,     NA,
  "TEST01", "PAT01", "PARAM02", 8.350, NA,       "Day 14",   NA,     "LOW",
  "TEST01", "PAT01", "PARAM03",    NA, "LOW",    "Baseline", "Y",    NA,
  "TEST01", "PAT01", "PARAM03",    NA, "LOW",    "Day 7",    NA,     NA,
  "TEST01", "PAT01", "PARAM03",    NA, "MEDIUM", "Day 14",   NA,     NA,
  "TEST01", "PAT01", "PARAM04",    NA, "HIGH",   "Baseline", "Y",    NA,
  "TEST01", "PAT01", "PARAM04",    NA, "HIGH",   "Day 7",    NA,     NA,
  "TEST01", "PAT01", "PARAM04",    NA, "MEDIUM", "Day 14",   NA,     NA
)

## Derive `BASE` variable from `AVAL`
derive_var_base(
  dataset,
  by_vars = exprs(USUBJID, PARAMCD),
  source_var = AVAL,
  new_var = BASE
)
#> # A tibble: 12 × 9
#>    STUDYID USUBJID PARAMCD  AVAL AVALC  AVISIT   ABLFL ANRIND  BASE
#>    <chr>   <chr>   <chr>   <dbl> <chr>  <chr>    <chr> <chr>  <dbl>
#>  1 TEST01  PAT01   PARAM01 10.1  NA     Baseline Y     NORMAL 10.1 
#>  2 TEST01  PAT01   PARAM01  9.7  NA     Day 7    NA    LOW    10.1 
#>  3 TEST01  PAT01   PARAM01 15.0  NA     Day 14   NA    HIGH   10.1 
#>  4 TEST01  PAT01   PARAM02  8.35 NA     Baseline Y     LOW     8.35
#>  5 TEST01  PAT01   PARAM02 NA    NA     Day 7    NA    NA      8.35
#>  6 TEST01  PAT01   PARAM02  8.35 NA     Day 14   NA    LOW     8.35
#>  7 TEST01  PAT01   PARAM03 NA    LOW    Baseline Y     NA     NA   
#>  8 TEST01  PAT01   PARAM03 NA    LOW    Day 7    NA    NA     NA   
#>  9 TEST01  PAT01   PARAM03 NA    MEDIUM Day 14   NA    NA     NA   
#> 10 TEST01  PAT01   PARAM04 NA    HIGH   Baseline Y     NA     NA   
#> 11 TEST01  PAT01   PARAM04 NA    HIGH   Day 7    NA    NA     NA   
#> 12 TEST01  PAT01   PARAM04 NA    MEDIUM Day 14   NA    NA     NA   

## Derive `BASEC` variable from `AVALC`
derive_var_base(
  dataset,
  by_vars = exprs(USUBJID, PARAMCD),
  source_var = AVALC,
  new_var = BASEC
)
#> # A tibble: 12 × 9
#>    STUDYID USUBJID PARAMCD  AVAL AVALC  AVISIT   ABLFL ANRIND BASEC
#>    <chr>   <chr>   <chr>   <dbl> <chr>  <chr>    <chr> <chr>  <chr>
#>  1 TEST01  PAT01   PARAM01 10.1  NA     Baseline Y     NORMAL NA   
#>  2 TEST01  PAT01   PARAM01  9.7  NA     Day 7    NA    LOW    NA   
#>  3 TEST01  PAT01   PARAM01 15.0  NA     Day 14   NA    HIGH   NA   
#>  4 TEST01  PAT01   PARAM02  8.35 NA     Baseline Y     LOW    NA   
#>  5 TEST01  PAT01   PARAM02 NA    NA     Day 7    NA    NA     NA   
#>  6 TEST01  PAT01   PARAM02  8.35 NA     Day 14   NA    LOW    NA   
#>  7 TEST01  PAT01   PARAM03 NA    LOW    Baseline Y     NA     LOW  
#>  8 TEST01  PAT01   PARAM03 NA    LOW    Day 7    NA    NA     LOW  
#>  9 TEST01  PAT01   PARAM03 NA    MEDIUM Day 14   NA    NA     LOW  
#> 10 TEST01  PAT01   PARAM04 NA    HIGH   Baseline Y     NA     HIGH 
#> 11 TEST01  PAT01   PARAM04 NA    HIGH   Day 7    NA    NA     HIGH 
#> 12 TEST01  PAT01   PARAM04 NA    MEDIUM Day 14   NA    NA     HIGH 

## Derive `BNRIND` variable from `ANRIND`
derive_var_base(
  dataset,
  by_vars = exprs(USUBJID, PARAMCD),
  source_var = ANRIND,
  new_var = BNRIND
)
#> # A tibble: 12 × 9
#>    STUDYID USUBJID PARAMCD  AVAL AVALC  AVISIT   ABLFL ANRIND BNRIND
#>    <chr>   <chr>   <chr>   <dbl> <chr>  <chr>    <chr> <chr>  <chr> 
#>  1 TEST01  PAT01   PARAM01 10.1  NA     Baseline Y     NORMAL NORMAL
#>  2 TEST01  PAT01   PARAM01  9.7  NA     Day 7    NA    LOW    NORMAL
#>  3 TEST01  PAT01   PARAM01 15.0  NA     Day 14   NA    HIGH   NORMAL
#>  4 TEST01  PAT01   PARAM02  8.35 NA     Baseline Y     LOW    LOW   
#>  5 TEST01  PAT01   PARAM02 NA    NA     Day 7    NA    NA     LOW   
#>  6 TEST01  PAT01   PARAM02  8.35 NA     Day 14   NA    LOW    LOW   
#>  7 TEST01  PAT01   PARAM03 NA    LOW    Baseline Y     NA     NA    
#>  8 TEST01  PAT01   PARAM03 NA    LOW    Day 7    NA    NA     NA    
#>  9 TEST01  PAT01   PARAM03 NA    MEDIUM Day 14   NA    NA     NA    
#> 10 TEST01  PAT01   PARAM04 NA    HIGH   Baseline Y     NA     NA    
#> 11 TEST01  PAT01   PARAM04 NA    HIGH   Day 7    NA    NA     NA    
#> 12 TEST01  PAT01   PARAM04 NA    MEDIUM Day 14   NA    NA     NA    
```
