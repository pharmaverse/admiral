# Derive First or Last Datetime from Multiple Sources

**\[deprecated\]**

The `derive_var_extreme_dtm()` function has been deprecated in favor of
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_vars_extreme_event.md).

Add the first or last datetime from multiple sources to the dataset,
e.g., the last known alive datetime (`LSTALVDTM`).

## Usage

``` r
derive_var_extreme_dtm(
  dataset,
  new_var,
  ...,
  source_datasets,
  mode,
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

- new_var:

  Name of variable to create

  Default value

  :   none

- ...:

  Source(s) of dates. One or more
  [`date_source()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/date_source.md)
  objects are expected.

  Default value

  :   none

- source_datasets:

  A named `list` containing datasets in which to search for the first or
  last date

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, the first date for each subject is
  selected. If `"last"` is specified, the last date for each subject is
  selected.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of expressions where the expressions are symbols as returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/reexport-exprs.md)
  is expected.

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

The input dataset with the new variable added.

## Details

The following steps are performed to create the output dataset:

1.  For each source dataset the observations as specified by the
    `filter` element are selected and observations where `date` is `NA`
    are removed. Then for each patient the first or last observation
    (with respect to `date` and `mode`) is selected.

2.  The new variable is set to the variable or expression specified by
    the `date` element. If this is a date variable (rather than
    datetime), then the time is imputed as `"00:00:00"`.

3.  The variables specified by the `set_values_to` element are added.

4.  The selected observations of all source datasets are combined into a
    single dataset.

5.  For each patient the first or last observation (with respect to the
    new variable and `mode`) from the single dataset is selected and the
    new variable is merged to the input dataset.

## See also

[`date_source()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/date_source.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_var_extreme_dt.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_vars_merged.md)

Other deprecated:
[`date_source()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_var_extreme_dt.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/3163-explore-use-of-dplyr-120/reference/get_summary_records.md)
