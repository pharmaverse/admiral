# Derive ATC Class Variables

Add Anatomical Therapeutic Chemical class variables from `FACM` to
`ADCM`.

**Note:** This is a wrapper function for the more generic
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_transposed.md).

## Usage

``` r
derive_vars_atc(
  dataset,
  dataset_facm,
  by_vars = exprs(!!!get_admiral_option("subject_keys"), CMREFID = FAREFID),
  id_vars = NULL,
  value_var = FASTRESC
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- dataset_facm:

  FACM dataset

  The variables specified by the `by_vars`, `id_vars`, and `value_var`
  arguments and `FATESTCD` are required. The variables `by_vars`,
  `id_vars`, and `FATESTCD` must be a unique key.

  Default value

  :   none

- by_vars:

  Grouping variables

  Keys used to merge `dataset_facm` with `dataset`.

  Default value

  :   `exprs(!!!get_admiral_option("subject_keys"), CMREFID = FAREFID)`

- id_vars:

  ID variables

  Variables (excluding by_vars) that uniquely identify each observation
  in `dataset_merge`.

  Default value

  :   `NULL`

- value_var:

  The variable of `dataset_facm` containing the values of the transposed
  variables

  Default value

  :   `FASTRESC`

## Value

The input dataset with ATC variables added

## See also

[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_transposed.md)

OCCDS Functions:
[`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_trtemfl.md),
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_query.md)

## Examples

``` r
library(tibble)

cm <- tribble(
  ~STUDYID,  ~USUBJID,       ~CMGRPID, ~CMREFID,  ~CMDECOD,
  "STUDY01", "BP40257-1001", "14",     "1192056", "PARACETAMOL",
  "STUDY01", "BP40257-1001", "18",     "2007001", "SOLUMEDROL",
  "STUDY01", "BP40257-1002", "19",     "2791596", "SPIRONOLACTONE"
)
facm <- tribble(
  ~STUDYID,  ~USUBJID,       ~FAGRPID, ~FAREFID,  ~FATESTCD,  ~FASTRESC,
  "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC1CD", "N",
  "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC2CD", "N02",
  "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC3CD", "N02B",
  "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC4CD", "N02BE",
  "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC1CD", "D",
  "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC2CD", "D10",
  "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC3CD", "D10A",
  "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC4CD", "D10AA",
  "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC1CD", "D",
  "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC2CD", "D07",
  "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC3CD", "D07A",
  "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC4CD", "D07AA",
  "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC1CD", "H",
  "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC2CD", "H02",
  "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC3CD", "H02A",
  "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC4CD", "H02AB",
  "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC1CD", "C",
  "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC2CD", "C03",
  "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC3CD", "C03D",
  "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC4CD", "C03DA"
)

derive_vars_atc(cm, facm, id_vars = exprs(FAGRPID))
#> # A tibble: 5 × 9
#>   STUDYID USUBJID      CMGRPID CMREFID CMDECOD       ATC1CD ATC2CD ATC3CD ATC4CD
#>   <chr>   <chr>        <chr>   <chr>   <chr>         <chr>  <chr>  <chr>  <chr> 
#> 1 STUDY01 BP40257-1001 14      1192056 PARACETAMOL   N      N02    N02B   N02BE 
#> 2 STUDY01 BP40257-1001 18      2007001 SOLUMEDROL    D      D10    D10A   D10AA 
#> 3 STUDY01 BP40257-1001 18      2007001 SOLUMEDROL    D      D07    D07A   D07AA 
#> 4 STUDY01 BP40257-1001 18      2007001 SOLUMEDROL    H      H02    H02A   H02AB 
#> 5 STUDY01 BP40257-1002 19      2791596 SPIRONOLACTO… C      C03    C03D   C03DA 
```
