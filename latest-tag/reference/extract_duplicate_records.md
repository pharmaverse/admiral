# Extract Duplicate Records

Extract Duplicate Records

## Usage

``` r
extract_duplicate_records(dataset, by_vars = NULL)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Defines groups of records in which to look for duplicates. If omitted,
  all variables in the input dataset are used in the by group.

  **Note:** Omitting `by_vars` will increase the function's run-time, so
  it is recommended to specify the necessary grouping variables for
  large datasets whenever possible.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

## Value

A `data.frame` of duplicate records within `dataset`

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/admiral-package.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/format.basket_select.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/signal_duplicate_records.md)

## Examples

``` r
data(admiral_adsl)

# Duplicate the first record
adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)

extract_duplicate_records(adsl, exprs(USUBJID))
#> # A tibble: 2 × 54
#>   USUBJID     STUDYID  SUBJID RFSTDTC RFENDTC RFXSTDTC RFXENDTC RFICDTC RFPENDTC
#>   <chr>       <chr>    <chr>  <chr>   <chr>   <chr>    <chr>    <chr>   <chr>   
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
