# Create Summary Records

**\[deprecated\]** The `get_summary_records()` has been deprecated in
favor of
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_summary_records.md)
(call it with the `dataset_add` argument and without the `dataset`
argument).

It is not uncommon to have an analysis need whereby one needs to derive
an analysis value (`AVAL`) from multiple records. The ADaM basic dataset
structure variable `DTYPE` is available to indicate when a new derived
records has been added to a dataset.

## Usage

``` r
get_summary_records(dataset, by_vars, filter = NULL, set_values_to = NULL)
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

  Variables to consider for generation of groupwise summary records.

  Default value

  :   none

- filter:

  Filter condition as logical expression to apply during summary
  calculation. By default, filtering expressions are computed within
  `by_vars` as this will help when an aggregating, lagging, or ranking
  function is involved.

  For example,

  - `filter_rows = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all
    AVAL values greater than mean of AVAL with in `by_vars`.

  - `filter_rows = (dplyr::n() > 2)` will filter n count of `by_vars`
    greater than 2.

  Default value

  :   `NULL`

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations.

  Set a list of variables to some specified value for the new records

  - LHS refer to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, an expression or NA. If summary
    functions are used, the values are summarized by the variables
    specified for `by_vars`.

  For example:

        set_values_to = exprs(
          AVAL = sum(AVAL),
          PARAMCD = "TDOSE",
          PARCAT1 = "OVERALL"
        )

  Default value

  :   `NULL`

## Value

A data frame of derived records.

## Details

This function only creates derived observations and does not append them
to the original dataset observations. If you would like to this instead,
see the
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_summary_records.md)
function.

## See also

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_summary_records.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_vars_merged_summary.md)

Other deprecated:
[`date_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/3170-adlb-template-paramcd-ck-labelled-creatinine-kinase-source-lbtest-is-creatine-kinase/reference/dthcaus_source.md)
