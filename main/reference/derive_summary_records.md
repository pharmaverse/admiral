# Add New Records Within By Groups Using Aggregation Functions

It is not uncommon to have an analysis need whereby one needs to derive
an analysis value (`AVAL`) from multiple records. The ADaM basic dataset
structure variable `DTYPE` is available to indicate when a new derived
records has been added to a dataset, if the derivation deviates from the
standard derivation of the parameter.

## Usage

``` r
derive_summary_records(
  dataset = NULL,
  dataset_add,
  dataset_ref = NULL,
  by_vars,
  filter_add = NULL,
  constant_values = NULL,
  set_values_to,
  missing_values = NULL
)
```

## Arguments

- dataset:

  Input dataset

  If the argument is not specified (or set to `NULL`), a new dataset is
  created. Otherwise, the new records are appended to the specified
  dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- dataset_add:

  Additional dataset

  The variables specified for `by_vars` are expected. Observations from
  the specified dataset are going to be used to calculate and added as
  new records to the input dataset (`dataset`).

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_ref:

  Reference dataset

  The variables specified for `by_vars` are expected. For each
  observation of the specified dataset a new observation is added to the
  input dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- by_vars:

  Grouping variables

  Variables to consider for generation of groupwise summary records.
  Providing the names of variables in
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  will create a groupwise summary and generate summary records for the
  specified groups.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- filter_add:

  Filter condition as logical expression to apply during summary
  calculation. By default, filtering expressions are computed within
  `by_vars` as this will help when an aggregating, lagging, or ranking
  function is involved.

  For example,

  - `filter_add = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all
    `AVAL` values greater than mean of `AVAL` with in `by_vars`.

  - `filter_add = (dplyr::n() > 2)` will filter n count of `by_vars`
    greater than 2.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- constant_values:

  Constant variables to set

  The specified variables are set to the specified values for all new
  summary records, including those with data in `dataset_add` and those
  with no data imputed using `dataset_ref` and `missing_values`.

  Set a list of variables to some specified value for the new records

  - LHS refer to a variable.

  - RHS refers to the values to set to the variable. This can be an
    expression.

  Permitted values

  :   list of named expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = "ml")`

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
    specified for `by_vars`. Any expression on the RHS must result in a
    single value per by group.

  For example:

        set_values_to = exprs(
          AVAL = sum(AVAL),
          DTYPE = "AVERAGE",
        )

  Permitted values

  :   list of named expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = "ml")`

  Default value

  :   none

- missing_values:

  Values for missing summary values

  For observations of the reference dataset (`dataset_ref`) which do not
  have a complete mapping defined by the summarization defined in
  `set_values_to`. Only variables specified for `set_values_to` can be
  specified for `missing_values`.

  Permitted values

  :   list of named expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = "ml")`

  Default value

  :   `NULL`

## Value

A data frame with derived records appended to original dataset.

## Details

For the newly derived records, only variables specified within `by_vars`
or `set_values_to` will be populated. All other variables will be set to
`NA`.

## See also

[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md)

## Examples

### Data setup

The following examples use the ECG dataset below as a basis.

    library(tibble, warn.conflicts = FALSE)
    library(dplyr, warn.conflicts = FALSE)

    adeg <- tribble(
      ~USUBJID,   ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL,
      "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,
      "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,
      "XYZ-1001", "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,
      "XYZ-1001", "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,
      "XYZ-1001", "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,
      "XYZ-1001", "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,
      "XYZ-1001", "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,
      "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,
      "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 200,
      "XYZ-1002", "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,
      "XYZ-1002", "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,
      "XYZ-1002", "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402
    ) %>%
      mutate(ADTM = convert_dtc_to_dtm(EGDTC))

### Summarize one or more variables using summary functions

A derived record is generated for each subject, containing the mean of
the triplicate ECG interval values (`AVAL`) and the latest measurement's
time (`ADTM`) by using summary functions within the `set_values_to`
argument.

    derive_summary_records(
      adeg,
      dataset_add = adeg,
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        ADTM = max(ADTM),
        DTYPE = "AVERAGE"
      )
    ) %>%
      arrange(USUBJID, AVISIT)
    #> # A tibble: 17 × 7
    #>    USUBJID  PARAM            AVISIT   EGDTC       AVAL ADTM                DTYPE
    #>    <chr>    <chr>            <chr>    <chr>      <dbl> <dttm>              <chr>
    #>  1 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  385  2016-02-24 07:50:00 <NA>
    #>  2 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-24 07:52:00 <NA>
    #>  3 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  396  2016-02-24 07:56:00 <NA>
    #>  4 XYZ-1001 QTcF Int. (msec) Baseline <NA>        393. 2016-02-24 07:56:00 AVER…
    #>  5 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  393  2016-03-08 09:48:00 <NA>
    #>  6 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  388  2016-03-08 09:51:00 <NA>
    #>  7 XYZ-1001 QTcF Int. (msec) Visit 2  <NA>        390. 2016-03-08 09:51:00 AVER…
    #>  8 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  394  2016-03-22 10:48:00 <NA>
    #>  9 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-22 10:51:00 <NA>
    #> 10 XYZ-1001 QTcF Int. (msec) Visit 3  <NA>        398  2016-03-22 10:51:00 AVER…
    #> 11 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-22 07:58:00 <NA>
    #> 12 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  200  2016-02-22 07:58:00 <NA>
    #> 13 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  392  2016-02-22 08:01:00 <NA>
    #> 14 XYZ-1002 QTcF Int. (msec) Baseline <NA>        330. 2016-02-22 08:01:00 AVER…
    #> 15 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  414  2016-03-24 10:53:00 <NA>
    #> 16 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-24 10:56:00 <NA>
    #> 17 XYZ-1002 QTcF Int. (msec) Visit 3  <NA>        408  2016-03-24 10:56:00 AVER…

Functions such as [`all()`](https://rdrr.io/r/base/all.html) and
[`any()`](https://rdrr.io/r/base/any.html) are also often useful when
creating summary records. For instance, the above example can be
extended to flag which derived records were affected by outliers. Note
that the outlier flag is created before `AVAL` is set for the summary
record. Otherwise, referencing `AVAL` later on would pick up the `AVAL`
from the summary record rather than the source records.

    derive_summary_records(
      adeg,
      dataset_add = adeg,
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(
        OUTLIEFL = if_else(any(AVAL >= 500 | AVAL <= 300), "Y", "N"),
        AVAL = mean(AVAL, na.rm = TRUE),
        ADTM = max(ADTM),
        DTYPE = "AVERAGE"
      )
    ) %>%
      arrange(USUBJID, AVISIT)
    #> # A tibble: 17 × 8
    #>    USUBJID  PARAM          AVISIT EGDTC  AVAL ADTM                OUTLIEFL DTYPE
    #>    <chr>    <chr>          <chr>  <chr> <dbl> <dttm>              <chr>    <chr>
    #>  1 XYZ-1001 QTcF Int. (ms… Basel… 2016…  385  2016-02-24 07:50:00 <NA>     <NA>
    #>  2 XYZ-1001 QTcF Int. (ms… Basel… 2016…  399  2016-02-24 07:52:00 <NA>     <NA>
    #>  3 XYZ-1001 QTcF Int. (ms… Basel… 2016…  396  2016-02-24 07:56:00 <NA>     <NA>
    #>  4 XYZ-1001 QTcF Int. (ms… Basel… <NA>   393. 2016-02-24 07:56:00 N        AVER…
    #>  5 XYZ-1001 QTcF Int. (ms… Visit… 2016…  393  2016-03-08 09:48:00 <NA>     <NA>
    #>  6 XYZ-1001 QTcF Int. (ms… Visit… 2016…  388  2016-03-08 09:51:00 <NA>     <NA>
    #>  7 XYZ-1001 QTcF Int. (ms… Visit… <NA>   390. 2016-03-08 09:51:00 N        AVER…
    #>  8 XYZ-1001 QTcF Int. (ms… Visit… 2016…  394  2016-03-22 10:48:00 <NA>     <NA>
    #>  9 XYZ-1001 QTcF Int. (ms… Visit… 2016…  402  2016-03-22 10:51:00 <NA>     <NA>
    #> 10 XYZ-1001 QTcF Int. (ms… Visit… <NA>   398  2016-03-22 10:51:00 N        AVER…
    #> 11 XYZ-1002 QTcF Int. (ms… Basel… 2016…  399  2016-02-22 07:58:00 <NA>     <NA>
    #> 12 XYZ-1002 QTcF Int. (ms… Basel… 2016…  200  2016-02-22 07:58:00 <NA>     <NA>
    #> 13 XYZ-1002 QTcF Int. (ms… Basel… 2016…  392  2016-02-22 08:01:00 <NA>     <NA>
    #> 14 XYZ-1002 QTcF Int. (ms… Basel… <NA>   330. 2016-02-22 08:01:00 Y        AVER…
    #> 15 XYZ-1002 QTcF Int. (ms… Visit… 2016…  414  2016-03-24 10:53:00 <NA>     <NA>
    #> 16 XYZ-1002 QTcF Int. (ms… Visit… 2016…  402  2016-03-24 10:56:00 <NA>     <NA>
    #> 17 XYZ-1002 QTcF Int. (ms… Visit… <NA>   408  2016-03-24 10:56:00 N        AVER…

### Restricting source records (`filter_add`)

The `filter_add` argument can be used to restrict the records that are
being summarized. For instance, the mean of the triplicates above can be
computed only for the baseline records by passing
`filter_add = AVISIT == "Baseline"`.

    derive_summary_records(
      adeg,
      dataset_add = adeg,
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      filter_add = AVISIT == "Baseline",
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        DTYPE = "AVERAGE"
      )
    ) %>%
      arrange(USUBJID, AVISIT)
    #> # A tibble: 14 × 7
    #>    USUBJID  PARAM            AVISIT   EGDTC       AVAL ADTM                DTYPE
    #>    <chr>    <chr>            <chr>    <chr>      <dbl> <dttm>              <chr>
    #>  1 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  385  2016-02-24 07:50:00 <NA>
    #>  2 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-24 07:52:00 <NA>
    #>  3 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  396  2016-02-24 07:56:00 <NA>
    #>  4 XYZ-1001 QTcF Int. (msec) Baseline <NA>        393. NA                  AVER…
    #>  5 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  393  2016-03-08 09:48:00 <NA>
    #>  6 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  388  2016-03-08 09:51:00 <NA>
    #>  7 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  394  2016-03-22 10:48:00 <NA>
    #>  8 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-22 10:51:00 <NA>
    #>  9 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-22 07:58:00 <NA>
    #> 10 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  200  2016-02-22 07:58:00 <NA>
    #> 11 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  392  2016-02-22 08:01:00 <NA>
    #> 12 XYZ-1002 QTcF Int. (msec) Baseline <NA>        330. NA                  AVER…
    #> 13 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  414  2016-03-24 10:53:00 <NA>
    #> 14 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-24 10:56:00 <NA> 

Summary functions can also be used within `filter_add` to filter based
on conditions applied to the whole of the by group specified in
`by_vars`. For instance, the mean of the triplicates can be computed
only for by groups which do indeed contain three records by passing
`filter_add = n() > 2`.

    derive_summary_records(
      adeg,
      dataset_add = adeg,
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      filter_add = n() > 2,
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        DTYPE = "AVERAGE"
      )
    ) %>%
      arrange(USUBJID, AVISIT)
    #> # A tibble: 14 × 7
    #>    USUBJID  PARAM            AVISIT   EGDTC       AVAL ADTM                DTYPE
    #>    <chr>    <chr>            <chr>    <chr>      <dbl> <dttm>              <chr>
    #>  1 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  385  2016-02-24 07:50:00 <NA>
    #>  2 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-24 07:52:00 <NA>
    #>  3 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  396  2016-02-24 07:56:00 <NA>
    #>  4 XYZ-1001 QTcF Int. (msec) Baseline <NA>        393. NA                  AVER…
    #>  5 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  393  2016-03-08 09:48:00 <NA>
    #>  6 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  388  2016-03-08 09:51:00 <NA>
    #>  7 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  394  2016-03-22 10:48:00 <NA>
    #>  8 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-22 10:51:00 <NA>
    #>  9 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-22 07:58:00 <NA>
    #> 10 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  200  2016-02-22 07:58:00 <NA>
    #> 11 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  392  2016-02-22 08:01:00 <NA>
    #> 12 XYZ-1002 QTcF Int. (msec) Baseline <NA>        330. NA                  AVER…
    #> 13 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  414  2016-03-24 10:53:00 <NA>
    #> 14 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-24 10:56:00 <NA> 

### Adding records for groups not in source (`dataset_ref` and `missing_values`)

Adding records for groups which are not in the source data can be
achieved by specifying a reference dataset in the `dataset_ref`
argument. For example, specifying the input dataset `adeg_allparamvis`
(containing an extra `"Visit 2"` for patient `1002`) ensures a summary
record is derived for that visit as well. For these records, the values
of the analysis variables to be populated should be specified within the
`missing_values` argument. Here, `DTYPE = "PHANTOM"` was chosen as
`AVAL` is set to missing.

    adeg_allparamvis <- tribble(
      ~USUBJID,   ~PARAM,             ~AVISIT,
      "XYZ-1001", "QTcF Int. (msec)", "Baseline",
      "XYZ-1001", "QTcF Int. (msec)", "Visit 2",
      "XYZ-1001", "QTcF Int. (msec)", "Visit 3",
      "XYZ-1002", "QTcF Int. (msec)", "Baseline",
      "XYZ-1002", "QTcF Int. (msec)", "Visit 2",
      "XYZ-1002", "QTcF Int. (msec)", "Visit 3"
    )

    derive_summary_records(
      adeg,
      dataset_add = adeg,
      dataset_ref = adeg_allparamvis,
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        ADTM = max(ADTM),
        DTYPE = "AVERAGE"
      ),
      missing_values = exprs(
        AVAL = NA,
        ADTM = NA,
        DTYPE = "PHANTOM"
      )
    ) %>%
      arrange(USUBJID, AVISIT)
    #> # A tibble: 18 × 7
    #>    USUBJID  PARAM            AVISIT   EGDTC       AVAL ADTM                DTYPE
    #>    <chr>    <chr>            <chr>    <chr>      <dbl> <dttm>              <chr>
    #>  1 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  385  2016-02-24 07:50:00 <NA>
    #>  2 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-24 07:52:00 <NA>
    #>  3 XYZ-1001 QTcF Int. (msec) Baseline 2016-02-2…  396  2016-02-24 07:56:00 <NA>
    #>  4 XYZ-1001 QTcF Int. (msec) Baseline <NA>        393. 2016-02-24 07:56:00 AVER…
    #>  5 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  393  2016-03-08 09:48:00 <NA>
    #>  6 XYZ-1001 QTcF Int. (msec) Visit 2  2016-03-0…  388  2016-03-08 09:51:00 <NA>
    #>  7 XYZ-1001 QTcF Int. (msec) Visit 2  <NA>        390. 2016-03-08 09:51:00 AVER…
    #>  8 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  394  2016-03-22 10:48:00 <NA>
    #>  9 XYZ-1001 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-22 10:51:00 <NA>
    #> 10 XYZ-1001 QTcF Int. (msec) Visit 3  <NA>        398  2016-03-22 10:51:00 AVER…
    #> 11 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  399  2016-02-22 07:58:00 <NA>
    #> 12 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  200  2016-02-22 07:58:00 <NA>
    #> 13 XYZ-1002 QTcF Int. (msec) Baseline 2016-02-2…  392  2016-02-22 08:01:00 <NA>
    #> 14 XYZ-1002 QTcF Int. (msec) Baseline <NA>        330. 2016-02-22 08:01:00 AVER…
    #> 15 XYZ-1002 QTcF Int. (msec) Visit 2  <NA>         NA  NA                  PHAN…
    #> 16 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  414  2016-03-24 10:53:00 <NA>
    #> 17 XYZ-1002 QTcF Int. (msec) Visit 3  2016-03-2…  402  2016-03-24 10:56:00 <NA>
    #> 18 XYZ-1002 QTcF Int. (msec) Visit 3  <NA>        408  2016-03-24 10:56:00 AVER…

### Add constant values to derived and missing summary records.

The `constant_values` argument allows you to assign fixed, common values
to all summary records generated by the function. This is particularly
useful when you need to populate new information for observations
derived from the `dataset_add` as well as for new records created for
subjects present in `dataset_ref` but missing in `dataset_add`.

For example, if `ADSL` contains two subjects ("1" and "2"), but `ADAE`
only has adverse event information for "Subject 1",
`derive_summary_records` will:

1.  Create a summary record for "Subject 1" based on `ADAE`.

2.  Identify "Subject 2" (from `ADSL`) as having no corresponding
    records in `ADAE` and create a new record for it.

The `constant_values` argument ensures that all these generated summary
records (for both Subject 1 and the newly created record for Subject 2)
receive the same `PARAMCD`, `PARAM`, and `PARCAT1` values. Additionally,
`missing_values` is used to specifically set `AVAL` to 0 for "Subject 2"
when no adverse events are found.

    library(tibble)

    adsl <- tibble(USUBJID = c("1", "2"))

    adae <- tribble(
     ~USUBJID, ~AEDECOD,
     "1",      "Illness",
     "1",      "Pain"
    )

    derive_summary_records(
      dataset_add = adae,
      dataset_ref = adsl,
      by_vars = exprs(USUBJID),
      constant_values = exprs(
        PARAMCD = "AECOUNT",
        PARAM = "Number of adverse events",
        PARCAT1 = "Adverse events"
      ),
      set_values_to = exprs(
        AVAL = n_distinct(AEDECOD),
        SRCDOM = "ADAE"
      ),
      missing_values = exprs(
        AVAL = 0
      )
    )
    #> # A tibble: 2 × 6
    #>   USUBJID  AVAL SRCDOM PARAMCD PARAM                    PARCAT1
    #>   <chr>   <dbl> <chr>  <chr>   <chr>                    <chr>
    #> 1 1           2 ADAE   AECOUNT Number of adverse events Adverse events
    #> 2 2           0 <NA>   AECOUNT Number of adverse events Adverse events
