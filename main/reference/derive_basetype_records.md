# Derive Basetype Variable

Baseline Type `BASETYPE` is needed when there is more than one
definition of baseline for a given Analysis Parameter `PARAM` in the
same dataset. For a given parameter, if Baseline Value `BASE` or `BASEC`
are derived and there is more than one definition of baseline, then
`BASETYPE` must be non-null on all records of any type for that
parameter where either `BASE` or `BASEC` are also non-null. Each value
of `BASETYPE` refers to a definition of baseline that characterizes the
value of `BASE` on that row. Please see section 4.2.1.6 of the ADaM
Implementation Guide, version 1.3 for further background.

## Usage

``` r
derive_basetype_records(dataset, basetypes)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `basetypes` argument are expected to be
  in the dataset.

  Default value

  :   none

- basetypes:

  A *named* list of expressions created using the
  [`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  function

  The names corresponds to the values of the newly created `BASETYPE`
  variables and the expressions are used to subset the input dataset.

  Default value

  :   none

## Value

The input dataset with variable `BASETYPE` added

## Details

Adds the `BASETYPE` variable to a dataset and duplicates records based
upon the provided conditions.

For each element of `basetypes` the input dataset is subset based upon
the provided expression and the `BASETYPE` variable is set to the name
of the expression. Then, all subsets are stacked. Records which do not
match any condition are kept and `BASETYPE` is set to `NA`.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

bds <- tribble(
  ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
  "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
  "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
  "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
  "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1,
  "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4,
  "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9,
  "P02",    "RUN-IN",       "PARAM01",     1,  12.1,
  "P02",    "DOUBLE-BLIND", "PARAM01",     2,  10.2,
  "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.8,
  "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4,
  "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8
)

bds_with_basetype <- derive_basetype_records(
  dataset = bds,
  basetypes = exprs(
    "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
    "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
    "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
  )
)


# Below print statement will print all 23 records in the data frame
# bds_with_basetype
print(bds_with_basetype, n = Inf)
#> # A tibble: 23 × 6
#>    USUBJID EPOCH        PARAMCD  ASEQ  AVAL BASETYPE    
#>    <chr>   <chr>        <chr>   <dbl> <dbl> <chr>       
#>  1 P01     RUN-IN       PARAM01     1  10   RUN-IN      
#>  2 P01     RUN-IN       PARAM01     2   9.8 RUN-IN      
#>  3 P01     DOUBLE-BLIND PARAM01     3   9.2 RUN-IN      
#>  4 P01     DOUBLE-BLIND PARAM01     4  10.1 RUN-IN      
#>  5 P01     OPEN-LABEL   PARAM01     5  10.4 RUN-IN      
#>  6 P01     OPEN-LABEL   PARAM01     6   9.9 RUN-IN      
#>  7 P02     RUN-IN       PARAM01     1  12.1 RUN-IN      
#>  8 P02     DOUBLE-BLIND PARAM01     2  10.2 RUN-IN      
#>  9 P02     DOUBLE-BLIND PARAM01     3  10.8 RUN-IN      
#> 10 P02     OPEN-LABEL   PARAM01     4  11.4 RUN-IN      
#> 11 P02     OPEN-LABEL   PARAM01     5  10.8 RUN-IN      
#> 12 P01     DOUBLE-BLIND PARAM01     3   9.2 DOUBLE-BLIND
#> 13 P01     DOUBLE-BLIND PARAM01     4  10.1 DOUBLE-BLIND
#> 14 P01     OPEN-LABEL   PARAM01     5  10.4 DOUBLE-BLIND
#> 15 P01     OPEN-LABEL   PARAM01     6   9.9 DOUBLE-BLIND
#> 16 P02     DOUBLE-BLIND PARAM01     2  10.2 DOUBLE-BLIND
#> 17 P02     DOUBLE-BLIND PARAM01     3  10.8 DOUBLE-BLIND
#> 18 P02     OPEN-LABEL   PARAM01     4  11.4 DOUBLE-BLIND
#> 19 P02     OPEN-LABEL   PARAM01     5  10.8 DOUBLE-BLIND
#> 20 P01     OPEN-LABEL   PARAM01     5  10.4 OPEN-LABEL  
#> 21 P01     OPEN-LABEL   PARAM01     6   9.9 OPEN-LABEL  
#> 22 P02     OPEN-LABEL   PARAM01     4  11.4 OPEN-LABEL  
#> 23 P02     OPEN-LABEL   PARAM01     5  10.8 OPEN-LABEL  

count(bds_with_basetype, BASETYPE, name = "Number of Records")
#> # A tibble: 3 × 2
#>   BASETYPE     `Number of Records`
#>   <chr>                      <int>
#> 1 DOUBLE-BLIND                   8
#> 2 OPEN-LABEL                     4
#> 3 RUN-IN                        11

# An example where all parameter records need to be included for 2 different
# baseline type derivations (such as LAST and WORST)
bds <- tribble(
  ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
  "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
  "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
  "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
  "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1
)

bds_with_basetype <- derive_basetype_records(
  dataset = bds,
  basetypes = exprs(
    "LAST" = TRUE,
    "WORST" = TRUE
  )
)

print(bds_with_basetype, n = Inf)
#> # A tibble: 8 × 6
#>   USUBJID EPOCH        PARAMCD  ASEQ  AVAL BASETYPE
#>   <chr>   <chr>        <chr>   <dbl> <dbl> <chr>   
#> 1 P01     RUN-IN       PARAM01     1  10   LAST    
#> 2 P01     RUN-IN       PARAM01     2   9.8 LAST    
#> 3 P01     DOUBLE-BLIND PARAM01     3   9.2 LAST    
#> 4 P01     DOUBLE-BLIND PARAM01     4  10.1 LAST    
#> 5 P01     RUN-IN       PARAM01     1  10   WORST   
#> 6 P01     RUN-IN       PARAM01     2   9.8 WORST   
#> 7 P01     DOUBLE-BLIND PARAM01     3   9.2 WORST   
#> 8 P01     DOUBLE-BLIND PARAM01     4  10.1 WORST   

count(bds_with_basetype, BASETYPE, name = "Number of Records")
#> # A tibble: 2 × 2
#>   BASETYPE `Number of Records`
#>   <chr>                  <int>
#> 1 LAST                       4
#> 2 WORST                      4
```
