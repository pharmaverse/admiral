# Summarize the Type-Specific Content of a Dataset

Collect the information
[`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)
reports for a dataset of class `ADSL`, `BDS`, `OCCDS`, or `TTE` (or
their medical device counterparts). An element is only returned when the
variables it needs are in the dataset. This matters because a class set
with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
is not checked against the variables, so e.g. a dataset declared as
`"TTE"` need not have `CNSR`.

## Usage

``` r
summarize_adsl(dataset, subject_vars)

summarize_bds(dataset, unit_vars, unit = "subjects")

summarize_occds(dataset, unit_vars)

summarize_tte(dataset)
```

## Arguments

- dataset:

  A data frame

  Default value

  :   none

- subject_vars:

  The variables subjects are counted by

  Default value

  :   none

- unit_vars:

  The variables the per-parameter counts and the checks are per:
  `subject_vars`, or for a medical device dataset the variables devices
  are identified by

  Default value

  :   `NULL`

- unit:

  The name of the per-parameter count column, `"subjects"` or
  `"devices"`

  Default value

  :   `NULL`

## Value

A named list:

- `summarize_adsl()`: `trt_var` and `trt` (subjects per treatment arm),
  `pops` (named count of subjects with `"Y"` per population flag),
  `flow` (subjects randomized, treated, and completed), `trtdur`
  (median, minimum, and maximum of `TRTDURD`), `eosstt` (subjects per
  end of study status), `n_deaths` and `dthcaus` (subjects who died, and
  per cause of death), `age` (median, minimum, and maximum of `AGE`),
  `sex` (subjects per `SEX`).

- `summarize_bds()`: `params` (per `PARAMCD`: the `PARAM` values, the
  records, subjects (or devices), and analysis visits, and the missing
  values, minimum, median, and maximum of a numeric `AVAL`), `visits`
  (the analysis visits, in `AVISITN` order if available), `dtype`
  (records per `DTYPE`), `baseline_by` and `multiple_baselines` (the
  variables a baseline is unique within, and the number of their
  combinations with more than one `ABLFL == "Y"` record).

- `summarize_occds()`: `terms` (named count of distinct values per
  coding variable), `trtem` (treatment-emergent and total records),
  `sev_var` and `sev` (the severity or toxicity grade variable, and the
  records per value), `ser_var` and `serious` (the seriousness variable,
  and the serious and total records), `occ_flag_by` (per occurrence flag
  checked, the variables it is unique within), `occ_flag_dups` (per
  occurrence flag checked, the number of combinations of those with more
  than one `"Y"`), `occ_flags_unchecked` (the occurrence flags whose
  level is not known).

- `summarize_tte()`: `events` (events and censored per `PARAMCD`).

## Details

The ADSL items other than the treatment arms and population flags are
taken from the first record of each subject, so that duplicated subject
records – which the record structure check reports – do not inflate
them. Every ADSL item counts subjects, so none is returned for a dataset
without any of the variables subjects are counted by.

The OCCDS severity is the first populated of `ASEV`, `--SEV`, `ATOXGR`,
and `--TOXGR`, and the seriousness the first populated `--SER`. See
[`check_occurrence_flags()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/check_occurrence_flags.md)
for the occurrence flag check.

A BDS dataset is expected to have at most one baseline record
(`ABLFL == "Y"`) per subject and parameter, and per `BASETYPE`, `ATPT`,
and `APERIOD` where these are in the dataset: a baseline per timepoint
or per period (e.g. in a crossover study) is correct, not a duplicate.
More than one silently corrupts `BASE`, `CHG`, and `PCHG`, and the
record structure check cannot see it, as the records differ in their
visit. Fewer than one is legitimate, as not every parameter has a
baseline. Records which belong to no subject (see
[`filter_unit_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/filter_unit_records.md))
are not checked.

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
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/signal_duplicate_records.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md),
[`summarize_vs_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_vs_adsl.md)
