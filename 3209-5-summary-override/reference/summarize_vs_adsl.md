# Compare a Dataset with its ADSL

Collect the information
[`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)
reports when its `adsl` argument is supplied.

## Usage

``` r
summarize_vs_adsl(dataset, adsl, subject_vars, dataset_class)
```

## Arguments

- dataset:

  A data frame

  Default value

  :   none

- adsl:

  The `adsl` argument of
  [`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)

  Default value

  :   none

- subject_vars:

  The variables subjects are counted and matched by

  Default value

  :   none

- dataset_class:

  The ADaM class of `dataset`

  Default value

  :   none

## Value

`NULL` (with a warning) if the comparison is not possible, otherwise a
named list with the elements `n_adsl` and `n_common` (the subjects in
`adsl`, and those of them with records), `orphans` (the subjects with
records which are not in `adsl`, as a data frame of `subject_vars`),
`n_orphans` (their number), and depending on the variables:

- `saffl`: the subjects with records within the safety population, and
  the subjects in it (`covered` and `total`),

- `arm_var` and `arm`: the treatment arm variable, and per arm the
  subjects with records and the subjects in it,

- `incidence` and `incidence_denom`: for an `OCCDS` dataset, the
  subjects with a treatment-emergent record and the subjects they are
  counted out of (`subjects` and `total`), and which subjects these are
  (`"in the safety population"` or `"in ADSL"`).

## Details

An `adsl` which cannot be compared with – one without `subject_vars`, or
with more than one record per subject – is an error, as any comparison
would be wrong. A `dataset` which cannot be compared, by contrast, only
gives a warning, as the rest of the summary is still valid.

Records of `dataset` with a missing subject variable are left out, as
they belong to no subject (see
[`filter_unit_records()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/filter_unit_records.md)).

The treatment arm is the first populated of `TRT01A`, `TRT01P`, and
`ARM` in `adsl`: the actual treatment ahead of the planned one, as the
subjects with records are most often the treated ones.

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
[`summarize_adsl()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_admiral_df_type.md),
[`summarize_bds_params()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_bds_params.md),
[`summarize_distribution()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_distribution.md),
[`summarize_params_by()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_params_by.md),
[`summarize_structure()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_structure.md),
[`summarize_subject_flow()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summarize_subject_flow.md)

## Examples

``` r
library(tibble)

adsl <- tribble(
  ~STUDYID, ~USUBJID, ~SAFFL, ~TRT01A,
  "PILOT01", "01-701-1015", "Y", "Placebo",
  "PILOT01", "01-701-1028", "Y", "Xanomeline High Dose",
  "PILOT01", "01-701-1034", "N", NA
)

adae <- tribble(
  ~STUDYID, ~USUBJID, ~AEDECOD, ~TRTEMFL,
  "PILOT01", "01-701-1015", "HEADACHE", "Y",
  "PILOT01", "01-701-1015", "NAUSEA", "N",
  "PILOT01", "01-701-9999", "RASH", "Y"
)

# `01-701-9999` is not in ADSL, so it is an orphan
str(admiral:::summarize_vs_adsl(
  adae,
  adsl,
  subject_vars = "USUBJID",
  dataset_class = "OCCDS"
))
#> List of 9
#>  $ n_adsl         : int 3
#>  $ n_common       : int 1
#>  $ n_orphans      : int 1
#>  $ orphans        : tibble [1 × 1] (S3: tbl_df/tbl/data.frame)
#>   ..$ USUBJID: chr "01-701-9999"
#>  $ saffl          : Named int [1:2] 1 2
#>   ..- attr(*, "names")= chr [1:2] "covered" "total"
#>  $ arm_var        : chr "TRT01A"
#>  $ arm            : tibble [3 × 3] (S3: tbl_df/tbl/data.frame)
#>   ..$ TRT01A  : chr [1:3] "Placebo" "Xanomeline High Dose" NA
#>   ..$ subjects: int [1:3] 1 0 0
#>   ..$ total   : int [1:3] 1 1 1
#>  $ incidence      : Named int [1:2] 1 2
#>   ..- attr(*, "names")= chr [1:2] "subjects" "total"
#>  $ incidence_denom: chr "in the safety population"
```
