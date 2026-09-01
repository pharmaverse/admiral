# Derive Death Cause

**\[deprecated\]** The `derive_var_dthcaus()` function has been
deprecated in favor of
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_vars_extreme_event.md).

Derive death cause (`DTHCAUS`) and add traceability variables if
required.

## Usage

``` r
derive_var_dthcaus(
  dataset,
  ...,
  source_datasets,
  subject_keys = get_admiral_option("subject_keys")
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `subject_keys` argument are expected to
  be in the dataset.

  Default value

  :   none

- ...:

  Objects of class "dthcaus_source" created by
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/dthcaus_source.md).

  Default value

  :   none

- source_datasets:

  A named `list` containing datasets in which to search for the death
  cause

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of expressions where the expressions are symbols as returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/reexport-exprs.md)
  is expected.

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

The input dataset with `DTHCAUS` variable added.

## Details

This function derives `DTHCAUS` along with the user-defined traceability
variables, if required. If a subject has death info from multiple
sources, the one from the source with the earliest death date will be
used. If dates are equivalent, the first source will be kept, so the
user should provide the inputs in the preferred order.

## See also

[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/dthcaus_source.md)

Other deprecated:
[`date_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_param_extreme_record.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/get_summary_records.md)
