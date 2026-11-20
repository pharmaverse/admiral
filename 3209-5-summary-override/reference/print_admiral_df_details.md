# Print the Type-Specific Part of a Summary

Print the lines of
[`print.summary_admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.summary_admiral_df.md)
for the class-specific items collected by
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md)
and
[`summarize_bds()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md).
Each line is only printed if its item is present.

## Usage

``` r
print_adsl_details(d)

print_bds_details(d)

print_occds_details(d)
```

## Arguments

- d:

  The `details` element of a `summary_admiral_df` object

  Default value

  :   none

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
[`print_summary_checks()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_checks.md),
[`print_summary_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_structure.md),
[`print_summary_table()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_summary_table.md),
[`print_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_vs_adsl.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/signal_duplicate_records.md),
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md),
[`summarize_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_vs_adsl.md)
