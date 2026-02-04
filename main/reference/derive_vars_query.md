# Derive Query Variables

Derive Query Variables

## Usage

``` r
derive_vars_query(dataset, dataset_queries)
```

## Arguments

- dataset:

  Input dataset

  Default value

  :   none

- dataset_queries:

  A dataset containing required columns `PREFIX`, `GRPNAME`, `SRCVAR`,
  `TERMCHAR` and/or `TERMNUM`, and optional columns `GRPID`, `SCOPE`,
  `SCOPEN`.

  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  can be used to create the dataset.

  Default value

  :   none

## Value

The input dataset with query variables derived.

## Details

This function can be used to derive CDISC variables such as `SMQzzNAM`,
`SMQzzCD`, `SMQzzSC`, `SMQzzSCN`, and `CQzzNAM` in ADAE and ADMH, and
variables such as `SDGzzNAM`, `SDGzzCD`, and `SDGzzSC` in ADCM. An
example usage of this function can be found in the
[`vignette("occds")`](https:/pharmaverse.github.io/admiral/main/articles/occds.md).

A query dataset is expected as an input to this function. See the
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/main/articles/queries_dataset.md)
for descriptions, or call `data("queries")` for an example of a query
dataset.

For each unique element in `PREFIX`, the corresponding "NAM" variable
will be created. For each unique `PREFIX`, if `GRPID` is not "" or NA,
then the corresponding "CD" variable is created; similarly, if `SCOPE`
is not "" or NA, then the corresponding "SC" variable will be created;
if `SCOPEN` is not "" or NA, then the corresponding "SCN" variable will
be created.

For each record in `dataset`, the "NAM" variable takes the value of
`GRPNAME` if the value of `TERMCHAR` or `TERMNUM` in `dataset_queries`
matches the value of the respective SRCVAR in `dataset`. Note that
`TERMCHAR` in `dataset_queries` dataset may be NA only when `TERMNUM` is
non-NA and vice versa. The matching is case insensitive. The "CD", "SC",
and "SCN" variables are derived accordingly based on `GRPID`, `SCOPE`,
and `SCOPEN` respectively, whenever not missing.

## See also

[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)

OCCDS Functions:
[`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtemfl.md),
[`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)

## Examples

``` r
library(tibble)
data("queries")
adae <- tribble(
  ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~AELLTCD,
  "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
  3, "Alanine aminotransferase abnormal", NA_character_, NA_integer_,
  "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
  5, "Basedow's disease", NA_character_, 1L,
  "03", "2020-06-07 23:59:59", "SOME TERM",
  2, "Some query", "Some term", NA_integer_,
  "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
  7, "Alveolar proteinosis", NA_character_, NA_integer_
)
derive_vars_query(adae, queries)
#> # A tibble: 4 × 24
#>   USUBJID ASTDTM     AETERM AESEQ AEDECOD AELLT AELLTCD SMQ02NAM SMQ02CD SMQ02SC
#>   <chr>   <chr>      <chr>  <dbl> <chr>   <chr>   <int> <chr>      <int> <chr>  
#> 1 01      2020-06-0… ALANI…     3 Alanin… NA         NA NA            NA NA     
#> 2 02      2020-06-0… BASED…     5 Basedo… NA          1 NA            NA NA     
#> 3 03      2020-06-0… SOME …     2 Some q… Some…      NA NA            NA NA     
#> 4 05      2020-06-0… ALVEO…     7 Alveol… NA         NA NA            NA NA     
#> # ℹ 14 more variables: SMQ02SCN <int>, SMQ03NAM <chr>, SMQ03CD <int>,
#> #   SMQ03SC <chr>, SMQ03SCN <int>, SMQ05NAM <chr>, SMQ05CD <int>,
#> #   SMQ05SC <chr>, SMQ05SCN <int>, CQ01NAM <chr>, CQ04NAM <chr>, CQ04CD <int>,
#> #   CQ06NAM <chr>, CQ06CD <int>
```
