# Get Query Variables

Create a table for the input dataset which binds the necessary rows for
a
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)
call with the relevant `SRCVAR`, `TERM_NAME_ID` and a temporary index if
it is necessary

**Note:** This function is the first step performed in
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)
requested by some users to be present independently from it.

## Usage

``` r
get_vars_query(dataset, dataset_queries)
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

  [`create_query_data()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_query_data.md)
  can be used to create the dataset.

  Default value

  :   none

## Value

The processed query dataset with `SRCVAR` and `TERM_NAME_ID` so that
that can be merged to the input dataset to execute the derivations
outlined by `dataset_queries`.

## Details

This function can be used to derive CDISC variables such as `SMQzzNAM`,
`SMQzzCD`, `SMQzzSC`, `SMQzzSCN`, and `CQzzNAM` in ADAE and ADMH, and
variables such as `SDGzzNAM`, `SDGzzCD`, and `SDGzzSC` in ADCM. An
example usage of this function can be found in the
[`vignette("occds")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/occds.md).

A query dataset is expected as an input to this function. See the
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/queries_dataset.md)
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

[`create_query_data()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_query_data.md)

Utilities used within Derivation functions:
[`extract_unit()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/extract_unit.md),
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_not_mapped.md)

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
get_vars_query(adae, queries)
#> # A tibble: 2 × 21
#>   USUBJID ASTDTM AETERM AESEQ SMQ02NAM SMQ02CD SMQ02SC SMQ02SCN SMQ03NAM SMQ03CD
#>   <chr>   <chr>  <chr>  <dbl> <chr>      <int> <chr>      <int> <chr>      <int>
#> 1 02      2020-… BASED…     5 NA            NA NA            NA NA            NA
#> 2 05      2020-… ALVEO…     7 NA            NA NA            NA NA            NA
#> # ℹ 11 more variables: SMQ03SC <chr>, SMQ03SCN <int>, SMQ05NAM <chr>,
#> #   SMQ05CD <int>, SMQ05SC <chr>, SMQ05SCN <int>, CQ01NAM <chr>, CQ04NAM <chr>,
#> #   CQ04CD <int>, CQ06NAM <chr>, CQ06CD <int>
```
