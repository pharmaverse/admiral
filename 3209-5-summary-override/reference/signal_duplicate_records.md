# Signal Duplicate Records

Signal Duplicate Records

## Usage

``` r
signal_duplicate_records(
  dataset,
  by_vars,
  msg = paste("Dataset contains duplicate records", "with respect to",
    "{.var {replace_values_by_names(by_vars)}}"),
  cnd_type = "error",
  class = NULL
)
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

  Defines groups of records in which to look for duplicates.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- msg:

  The condition message

  Default value

  :   `paste("Dataset contains duplicate records", "with respect to", "{.var {replace_values_by_names(by_vars)}}")`

- cnd_type:

  Type of condition to signal when detecting duplicate records.

  Permitted values

  :   `"message"`, `"warning"`, or `"error"`

  Default value

  :   `"error"`

- class:

  Class of the condition

  The specified classes are added to the classes of the condition.
  `c("duplicate_records", "assert-admiral")` is always added.

  Default value

  :   `NULL`

## Value

No return value, called for side effects

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral-package.md),
[`check_occurrence_flags()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/check_occurrence_flags.md),
[`count_incidence()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_incidence.md),
[`count_multiple_baselines()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_multiple_baselines.md),
[`count_units()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_units.md),
[`count_values()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/count_values.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/extract_duplicate_records.md),
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
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md),
[`summarize_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_vs_adsl.md)

## Examples

``` r
data(admiral_adsl)

# Duplicate the first record
adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)

signal_duplicate_records(adsl, exprs(USUBJID), cnd_type = "message")
#> Dataset contains duplicate records with respect to `USUBJID`
#> ℹ Run `admiral::get_duplicates_dataset()` to access the duplicate records
```
