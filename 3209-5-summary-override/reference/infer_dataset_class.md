# Infer the ADaM Dataset Class

Attempts to classify a data frame as one of the common ADaM dataset
classes from its variables and records. It ignores any class declared
with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md):
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_class.md)
returns the declared class if there is one, and only falls back to this
function otherwise.

Note that this is the *class* of the dataset (`ADSL`, `BDS`, ...), not
its record structure: the latter is the set of key variables which
uniquely identify its records. Unlike the class, the keys are never
inferred; they are only ever declared with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).

## Usage

``` r
infer_dataset_class(dataset)
```

## Arguments

- dataset:

  A data frame

  Default value

  :   none

## Value

A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`,
`"ADDL"`, `"MDBDS"`, `"MDTTE"`, `"MDOCCDS"`, or `"other"`; or `NA` if
the dataset has no records and its variables do not decide the class.

## Details

A dataset is first classed as device data or subject data:

- **Device data** contains `SPDEVID`, the device identifier, populated
  (not `NA`) for every record. Its class is one of the medical device
  classes `"ADDL"`, `"MDBDS"`, `"MDTTE"`, and `"MDOCCDS"`.

- **Subject data** is otherwise data which contains all the subject keys
  (`get_admiral_option("subject_keys")`, by default `STUDYID` and
  `USUBJID`). Its class is one of `"ADSL"`, `"BDS"`, `"OCCDS"`, and
  `"TTE"`.

A dataset which is neither is `"other"`. This keeps a data frame which
merely shares variable names with an ADaM dataset – a coding dictionary
with `AETERM`/`AEDECOD` columns, a summary table with `PARAMCD` and
`AVAL` – from being taken for one.

The subject keys do not decide the class: in the ADaMIG-MD `STUDYID` is
required in a device dataset and `USUBJID` is conditional, so a device
dataset may identify the subject of every record. Nor does the presence
of `SPDEVID`, as subject datasets record it too, e.g. the device
involved in an adverse event; but there it is populated only for the
records which involve a device. A subject dataset in which every record
involves a device (e.g. an `ADAE` restricted to device-related events)
is therefore inferred as device data, so declare its class with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
instead.

Within each class the following precedence is used. The order matters,
because the classes are not mutually exclusive by variable content – a
`TTE` dataset contains `PARAMCD` just as a `BDS` one does, and an
occurrence dataset which has had a parameter added is structurally both.

1.  `"TTE"`/`"MDTTE"` – contains `PARAMCD`, `CNSR`, and `STARTDT`.
    Tested first because a time-to-event dataset is a `BDS` dataset by
    variable content; only the censoring and origin-date variables
    distinguish it.

2.  `"BDS"`/`"MDBDS"` – contains `PARAMCD` and `AVAL` or `AVALC`.

3.  `"OCCDS"`/`"MDOCCDS"` – does not contain `PARAMCD`, and contains a
    `--DECOD`/`--TERM` variable with a two-letter SDTM domain prefix
    (`AEDECOD`, `DETERM`, ...; a variable merely ending in `TERM` does
    not count), an `AOCCxxFL` occurrence flag, or `TRTEMFL`.

4.  `"ADSL"`/`"ADDL"` – does not contain `PARAMCD` and has one record
    per subject (per combination of the subject keys) or per device (per
    combination of `STUDYID` and `SPDEVID`). This is decided by the
    records: treatment variables are no evidence of ADSL, as
    `BDS`/`OCCDS` derivations merge them on from ADSL.

5.  `"other"` – none of the above, if there is at least one record.

6.  `NA` – none of the above.

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
