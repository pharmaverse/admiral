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
      [`exprs()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

## Value

A `data.frame` of duplicate records within `dataset`

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral-package.md),
[`check_occurrence_flags()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/check_occurrence_flags.md),
[`count_incidence()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_incidence.md),
[`count_multiple_baselines()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_multiple_baselines.md),
[`count_units()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_units.md),
[`count_values()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_values.md),
[`filter_unit_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/filter_unit_records.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/format.basket_select.md),
[`format_counts()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/format_counts.md),
[`format_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/format_distribution.md),
[`format_number()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/format_number.md),
[`format_percent()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/format_percent.md),
[`infer_class_by_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/infer_class_by_structure.md),
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/infer_dataset_class.md),
[`new_admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/new_admiral_df.md),
[`print_adsl_details()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_admiral_df_details.md),
[`print_summary_checks()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_checks.md),
[`print_summary_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_structure.md),
[`print_summary_table()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_table.md),
[`print_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_vs_adsl.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/signal_duplicate_records.md),
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_params_by()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_params_by.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md),
[`summarize_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_vs_adsl.md)

## Examples

``` r
data(admiral_adsl)

# Duplicate the first record
adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)

extract_duplicate_records(adsl, exprs(USUBJID))
#> # An admiral dataset (name: ADSL, class: ADSL, keys: STUDYID, USUBJID)
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
