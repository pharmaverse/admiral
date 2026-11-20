# Check the Occurrence Flags of an OCCDS Dataset

Counts, for each occurrence flag (`AOCC...FL`), the groups with more
than one `"Y"` record. An occurrence flag marks the first occurrence per
subject, or per subject and body system or dictionary term, so more than
one `"Y"` per group counts a subject twice in every "n (%) of subjects"
table built from it – which the record structure check cannot see, as
the records differ in other variables. In a medical device dataset, the
flags are per device instead.

## Usage

``` r
check_occurrence_flags(dataset, unit_vars)
```

## Arguments

- dataset:

  A data frame

  Default value

  :   none

- unit_vars:

  The variables subjects, or devices, are identified by

  Default value

  :   none

## Value

A named list with the elements `occ_flag_by`, `occ_flag_dups`, and
`occ_flags_unchecked`, see
[`summarize_occds()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md);
each only if not empty.

## Details

The group a flag is unique within is taken from its name, following the
ADaMIG for OCCDS: `AOCCFL` per subject, `AOCCSFL` per subject and body
system (`--BODSYS`), and `AOCCPFL` per subject and dictionary term
(`--DECOD`). The maximum intensity versions (`AOCCIFL`, `AOCCSIFL`,
`AOCCPIFL`) are unique within the same groups. `APERIOD` is added to
every group if it is populated, as multi-period designs (e.g. vaccine
studies) flag the first occurrence per period.

Any other flag, e.g. a sponsor-defined `AOCC01FL`, is not checked, as
its group cannot be told from its name; neither is a flag whose body
system or dictionary term variable is not in the dataset. Fewer than one
`"Y"` per group is not a defect, as the flags are often restricted to a
subset of the records, such as the treatment-emergent ones. Records
which belong to no subject or device (see
[`filter_unit_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/filter_unit_records.md))
are not checked.

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral-package.md),
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
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/signal_duplicate_records.md),
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md),
[`summarize_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_vs_adsl.md)
