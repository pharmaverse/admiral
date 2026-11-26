# Construct an `admiral_df` Object

Low-level constructor of the `admiral_df` class: adds the class to a
data frame (unless it is already present), preserving the existing
classes such as `tbl_df`/`data.frame`, and stores any metadata supplied.
The arguments are not checked, so this is what admiral's `derive_*()`
and related dataset-producing functions apply to their output. Users
call
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md),
which validates its arguments before calling this.

## Usage

``` r
new_admiral_df(
  x,
  dataset_keys = NULL,
  dataset_name = NULL,
  dataset_class = NULL
)
```

## Arguments

- x:

  A data frame (or `NULL`)

  Default value

  :   none

- dataset_keys:

  Character vector of key variables, stored in the
  `"admiral_dataset_keys"` attribute

  Default value

  :   `NULL`

- dataset_name:

  Character scalar, stored in the `"admiral_dataset_name"` attribute

  Default value

  :   `NULL`

- dataset_class:

  Character scalar, stored in the `"admiral_dataset_class"` attribute

  Default value

  :   `NULL`

## Value

If `x` is `NULL`, `NULL` is returned unchanged. Otherwise `x` with the
supplied attributes set and `"admiral_df"` prepended to its class
attribute.

## Details

An argument left as `NULL` leaves the corresponding attribute untouched,
so that tagging a dataset which already carries metadata does not remove
it.

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
