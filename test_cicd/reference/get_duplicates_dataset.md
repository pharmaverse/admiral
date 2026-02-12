# Get Duplicate Records that Led to a Prior Error

Get Duplicate Records that Led to a Prior Error

## Usage

``` r
get_duplicates_dataset()
```

## Value

A `data.frame` or `NULL`

## Details

Many `{admiral}` function check that the input dataset contains only one
record per `by_vars` group and throw an error otherwise. The
`get_duplicates_dataset()` function allows one to retrieve the duplicate
records that lead to an error.

Note that the function always returns the dataset of duplicates from the
last error that has been thrown in the current R session. Thus, after
restarting the R sessions `get_duplicates_dataset()` will return `NULL`
and after a second error has been thrown, the dataset of the first error
can no longer be accessed (unless it has been saved in a variable).

## See also

Utilities for Dataset Checking:
[`get_many_to_one_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_many_to_one_dataset.md),
[`get_one_to_many_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_one_to_many_dataset.md)

## Examples

``` r
data(admiral_adsl)

# Duplicate the first record
adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)

signal_duplicate_records(adsl, exprs(USUBJID), cnd_type = "warning")
#> Warning: Dataset contains duplicate records with respect to `USUBJID`
#> ℹ Run `admiral::get_duplicates_dataset()` to access the duplicate records

get_duplicates_dataset()
#> Duplicate records with respect to `USUBJID`.
#> # A tibble: 2 × 54
#>   USUBJID     STUDYID  SUBJID RFSTDTC RFENDTC RFXSTDTC RFXENDTC RFICDTC RFPENDTC
#> * <chr>       <chr>    <chr>  <chr>   <chr>   <chr>    <chr>    <chr>   <chr>   
#> 1 01-701-1015 CDISCPI… 1015   2014-0… 2014-0… 2014-01… 2014-07… NA      2014-07…
#> 2 01-701-1015 CDISCPI… 1015   2014-0… 2014-0… 2014-01… 2014-07… NA      2014-07…
#> # ℹ 45 more variables: DTHDTC <chr>, DTHFL <chr>, SITEID <chr>, AGE <dbl>,
#> #   AGEU <chr>, SEX <chr>, RACE <chr>, ETHNIC <chr>, ARMCD <chr>, ARM <chr>,
#> #   ACTARMCD <chr>, ACTARM <chr>, COUNTRY <chr>, DMDTC <chr>, DMDY <dbl>,
#> #   TRT01P <chr>, TRT01A <chr>, TRTSDTM <dttm>, TRTSTMF <chr>, TRTEDTM <dttm>,
#> #   TRTETMF <chr>, TRTSDT <date>, TRTEDT <date>, TRTDURD <dbl>, SCRFDT <date>,
#> #   EOSDT <date>, EOSSTT <chr>, FRVDT <date>, RANDDT <date>, DTHDT <date>,
#> #   DTHDTF <chr>, DTHADY <dbl>, LDDTHELD <dbl>, DTHCAUS <chr>, DTHDOM <chr>, …
```
