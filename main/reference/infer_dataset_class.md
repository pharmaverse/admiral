# Infer the ADaM Dataset Class

Attempts to classify a data frame as one of the common ADaM dataset
classes from its variables and records. It ignores any class declared
with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md):
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/get_dataset_class.md)
returns the declared class if there is one, and only falls back to this
function otherwise.

Note that this is the *class* of the dataset (`ADSL`, `BDS`, ...), not
its record structure: the latter is the set of key variables which
uniquely identify its records. Unlike the class, the keys are never
inferred; they are only ever declared with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md).

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
  (neither `NA` nor `""`) for every record. Its class is one of the
  medical device classes `"ADDL"`, `"MDBDS"`, `"MDTTE"`, and
  `"MDOCCDS"`.

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
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md)
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
[`admiral-package`](https:/pharmaverse.github.io/admiral/main/reference/admiral-package.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/extract_duplicate_records.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/format.basket_select.md),
[`infer_class_by_structure()`](https:/pharmaverse.github.io/admiral/main/reference/infer_class_by_structure.md),
[`new_admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/new_admiral_df.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/signal_duplicate_records.md)
