# Add the First or Last Observation for Each By Group as New Records

Add the first or last observation for each by group as new observations.
The new observations can be selected from the additional dataset. This
function can be used for adding the maximum or minimum value as a
separate visit. All variables of the selected observation are kept. This
distinguishes `derive_extreme_records()` from
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md),
where only the by variables are populated for the new records.

## Usage

``` r
derive_extreme_records(
  dataset = NULL,
  dataset_add,
  dataset_ref = NULL,
  by_vars = NULL,
  order = NULL,
  mode = NULL,
  filter_add = NULL,
  check_type = "warning",
  exist_flag = NULL,
  true_value = "Y",
  false_value = NA_character_,
  keep_source_vars = exprs(everything()),
  set_values_to
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

  The additional dataset, which determines the by groups returned in the
  input dataset, based on the groups that exist in this dataset after
  being subset by `filter_add`.

  The variables specified in the `by_vars` and `filter_add` parameters
  are expected in this dataset. If `mode` and `order` are specified, the
  first or last observation within each by group, defined by `by_vars`,
  is selected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_ref:

  Reference dataset

  The variables specified for `by_vars` are expected. For each
  observation of the specified dataset a new observation is added to the
  input dataset.

  For records which are added from `dataset_ref` because there are no
  records in `dataset_add` for the by group only those variables are
  kept which are also in `dataset_add` (and are included in
  `keep_source_vars`).

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- by_vars:

  Grouping variables

  If `dataset_ref` is specified, this argument must be specified.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- order:

  Sort order

  Within each by group the observations are ordered by the specified
  order.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, the first observation of each by group is
  added to the input dataset. If `"last"` is specified, the last
  observation of each by group is added to the input dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   `NULL`

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations in `dataset_add` fulfilling the specified condition
  are considered.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the (restricted) additional dataset are
  not unique with respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

- exist_flag:

  Existence flag

  The specified variable is added to the output dataset.

  For by groups with at least one observation in the additional dataset
  (`dataset_add`) `exist_flag` is set to the value specified by the
  `true_value` argument.

  For all other by groups `exist_flag` is set to the value specified by
  the `false_value` argument.

  Permitted values

  :   Variable name

  Default value

  :   `NULL`

- true_value:

  True value

  For new observations selected from the additional dataset
  (`dataset_add`), `exist_flag` is set to the specified value.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"Y"`

- false_value:

  False value

  For new observations not selected from the additional dataset
  (`dataset_add`), `exist_flag` is set to the specified value.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- keep_source_vars:

  Variables to be kept in the new records

  A named list or tidyselect expressions created by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining the variables to be kept for the new records. The variables
  specified for `by_vars` and `set_values_to` need not be specified here
  as they are kept automatically.

  Permitted values

  :   list of variables or tidyselect expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(DTHDT, starts_with("AST"))` or `exprs(everything)`

  Default value

  :   `exprs(everything())`

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
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = "ml")`

  Default value

  :   none

## Value

The input dataset with the first or last observation of each by group
added as new observations.

## Details

1.  The additional dataset (`dataset_add`) is restricted as specified by
    the `filter_add` argument.

2.  For each group (with respect to the variables specified for the
    `by_vars` argument) the first or last observation (with respect to
    the order specified for the `order` argument and the mode specified
    for the `mode` argument) is selected.

3.  If `dataset_ref` is specified, observations which are in
    `dataset_ref` but not in the selected records are added. Variables
    that are common across `dataset_ref`, `dataset_add` and
    `keep_source_vars()` are also populated for the new observations.

4.  The variables specified by the `set_values_to` argument are added to
    the selected observations.

5.  The variables specified by the `keep_source_vars` argument are
    selected along with the variables specified in `by_vars` and
    `set_values_to` arguments.

6.  The observations are added to input dataset (`dataset`). If no input
    dataset is provided, a new dataset is created.

## See also

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

## Examples

### Add last/first record as new record

For each subject the last record should be added as a new visit.

- The source dataset for the new records is specified by the
  `dataset_add` argument. Here it is the same as the input dataset.

- The groups for which new records are added are specified by the
  `by_vars` argument. Here for each *subject* a record should be added.
  Thus `by_vars = exprs(USUBJID)` is specified.

- As there are multiple records per subject, the `mode` and `order`
  arguments are specified to request that the *last* record is selected
  when the records are sorted by visit (`AVISITN`). The records are
  sorted by each by group (`by_vars`) separately, i.e., it is not
  necessary to add the variables from `by_vars` to `order`.

- To avoid duplicates in the output dataset the `set_values_to` argument
  is specified to set the visit (`AVISIT`) to a special value for the
  *new* records.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)
    library(lubridate, warn.conflicts = FALSE)

    adlb <- tribble(
      ~USUBJID, ~AVISITN, ~AVAL,
      "1",      1,          113,
      "1",      2,          111,
      "2",      1,          101,
      "2",      2,           NA,
      "3",      1,           NA,
    )

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      by_vars = exprs(USUBJID),
      order = exprs(AVISITN),
      mode = "last",
      set_values_to = exprs(
        AVISITN = 99
      )
    )
    #> # A tibble: 8 × 3
    #>   USUBJID AVISITN  AVAL
    #>   <chr>     <dbl> <dbl>
    #> 1 1             1   113
    #> 2 1             2   111
    #> 3 2             1   101
    #> 4 2             2    NA
    #> 5 3             1    NA
    #> 6 1            99   111
    #> 7 2            99    NA
    #> 8 3            99    NA

### Restricting source records (`filter_add`)

The source records can be restricted by the `filter_add` argument, e.g.,
to exclude visits with missing analysis value from selecting for the
last visit record:

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL),
      by_vars = exprs(USUBJID),
      order = exprs(AVISITN),
      mode = "last",
      set_values_to = exprs(
        AVISITN = 99
      )
    )
    #> # A tibble: 7 × 3
    #>   USUBJID AVISITN  AVAL
    #>   <chr>     <dbl> <dbl>
    #> 1 1             1   113
    #> 2 1             2   111
    #> 3 2             1   101
    #> 4 2             2    NA
    #> 5 3             1    NA
    #> 6 1            99   111
    #> 7 2            99   101

Please note that new records are added only for subjects in the
*restricted* source data. Therefore no new record is added for subject
`"3"`.

### Adding records for groups not in source (`dataset_ref`)

Adding records for groups which are not in the source data can be
achieved by specifying a reference dataset by the `dataset_ref`
argument. For example, specifying the input dataset for `dataset_ref`
below ensures that new records are added also for subject without a
valid analysis value:

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL),
      dataset_ref = adlb,
      by_vars = exprs(USUBJID),
      order = exprs(AVISITN),
      mode = "last",
      set_values_to = exprs(
        AVISITN = 99
      )
    )
    #> # A tibble: 8 × 3
    #>   USUBJID AVISITN  AVAL
    #>   <chr>     <dbl> <dbl>
    #> 1 1             1   113
    #> 2 1             2   111
    #> 3 2             1   101
    #> 4 2             2    NA
    #> 5 3             1    NA
    #> 6 1            99   111
    #> 7 2            99   101
    #> 8 3            99    NA

### Selecting variables for new records (`keep_source_vars`)

Which variables from the source dataset are kept for the new records can
be specified by the `keep_source_vars` argument. Variables specified by
`by_vars` or `set_values_to` don't need to be added to
`keep_source_vars` as these are always kept.

    adlb <- tribble(
      ~USUBJID, ~AVISIT,  ~AVAL, ~LBSEQ,
      "1",      "WEEK 1",   123,      1,
      "1",      "WEEK 2",   101,      2,
      "2",      "WEEK 1",    99,      1,
      "2",      "WEEK 2",   110,      2,
      "2",      "WEEK 3",    93,      3
    )

    derive_extreme_records(
      dataset_add = adlb,
      filter_add = !is.na(AVAL),
      by_vars = exprs(USUBJID),
      order = exprs(AVAL),
      mode = "first",
      keep_source_vars = exprs(AVAL),
      set_values_to = exprs(
        AVISIT = "MINIMUM"
      )
    )
    #> # A tibble: 2 × 3
    #>   USUBJID AVISIT   AVAL
    #>   <chr>   <chr>   <dbl>
    #> 1 1       MINIMUM   101
    #> 2 2       MINIMUM    93

### Handling duplicates (`check_type`)

The source records are checked regarding duplicates with respect to
`by_vars` and `order`. By default, a warning is issued if any duplicates
are found.

    adlb <- tribble(
      ~USUBJID, ~AVISIT,  ~AVAL,
      "1",      "WEEK 1",   123,
      "1",      "WEEK 2",   123,
      "2",      "WEEK 1",    99,
      "2",      "WEEK 2",   110,
      "2",      "WEEK 3",    93,
    )

    derive_extreme_records(
      dataset_add = adlb,
      filter_add = !is.na(AVAL),
      by_vars = exprs(USUBJID),
      order = exprs(AVAL),
      mode = "first",
      set_values_to = exprs(
        AVISIT = "MINIMUM"
      )
    )
    #> # A tibble: 2 × 3
    #>   USUBJID AVISIT   AVAL
    #>   <chr>   <chr>   <dbl>
    #> 1 1       MINIMUM   123
    #> 2 2       MINIMUM    93
    #> Warning: Dataset contains duplicate records with respect to `USUBJID` and `AVAL`
    #> i Run `admiral::get_duplicates_dataset()` to access the duplicate records

For investigating the issue, the dataset of the duplicate source records
can be obtained by calling
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_duplicates_dataset.md):

    get_duplicates_dataset()
    #> Duplicate records with respect to `USUBJID` and `AVAL`.
    #> # A tibble: 2 × 3
    #>   USUBJID  AVAL AVISIT
    #> * <chr>   <dbl> <chr>
    #> 1 1         123 WEEK 1
    #> 2 1         123 WEEK 2

Common options to solve the issue are:

- Restricting the source records by specifying/updating the `filter_add`
  argument.

- Specifying additional variables for `order`.

- Setting `check_type = "none"` to ignore any duplicates.

In this example it doesn't matter which of the records with the minimum
value is chosen because it doesn't affect the output dataset. Thus the
third option is used:

    derive_extreme_records(
      dataset_add = adlb,
      filter_add = !is.na(AVAL),
      by_vars = exprs(USUBJID),
      order = exprs(AVAL),
      mode = "first",
      check_type = "none",
      set_values_to = exprs(
        AVISIT = "MINIMUM"
      )
    )
    #> # A tibble: 2 × 3
    #>   USUBJID AVISIT   AVAL
    #>   <chr>   <chr>   <dbl>
    #> 1 1       MINIMUM   123
    #> 2 2       MINIMUM    93

### Flagging existence of source records (`exist_flag`, `true_value`, `false_value`)

If the existence of a source record should be flagged, the `exist_flag`
argument can be specified. The specified variable is set to `true_value`
if a source record exists. Otherwise, it is set to `false_value`.

The `dataset_ref` argument should be specified as otherwise *all* new
records originate from `dataset_add`, i.e., `exist_flag` would be set to
`true_value` for all records.

    adsl <- tribble(
      ~USUBJID, ~DTHDT,
      "1",      ymd("2022-05-13"),
      "2",      ymd(""),
      "3",      ymd("")
    )

    derive_extreme_records(
      dataset_ref = adsl,
      dataset_add = adsl,
      by_vars = exprs(USUBJID),
      filter_add = !is.na(DTHDT),
      exist_flag = AVALC,
      true_value = "Y",
      false_value = "N",
      set_values_to = exprs(
        PARAMCD = "DEATH",
        ADT = DTHDT
      )
    )
    #> # A tibble: 3 × 5
    #>   USUBJID PARAMCD ADT        DTHDT      AVALC
    #>   <chr>   <chr>   <date>     <date>     <chr>
    #> 1 1       DEATH   2022-05-13 2022-05-13 Y
    #> 2 2       DEATH   NA         NA         N
    #> 3 3       DEATH   NA         NA         N    

### Derive `DTYPE = "LOV"`

For each subject and parameter the last valid assessment (with respect
to `AVISITN` and `LBSEQ`) should be selected and added as a new record
to the input dataset. For the new records set `AVISIT = "PBL LAST"`,
`AVISITN = 99`, and `DTYPE = "LOV"`.

    adlb <- tribble(
      ~USUBJID, ~AVISIT,    ~AVISITN, ~PARAMCD, ~AVAL, ~LBSEQ,
      "1",      "BASELINE",        1, "ABC",      120,      1,
      "1",      "WEEK 1",          2, "ABC",      113,      2,
      "1",      "WEEK 1",          2, "ABC",      117,      3,
      "2",      "BASELINE",        1, "ABC",      101,      1,
      "2",      "WEEK 1",          2, "ABC",      101,      2,
      "2",      "WEEK 2",          3, "ABC",       95,      3,
      "1",      "BASELINE",        1, "DEF",       17,      1,
      "1",      "WEEK 1",          2, "DEF",       NA,      2,
      "1",      "WEEK 1",          2, "DEF",       13,      3,
      "2",      "BASELINE",        1, "DEF",        9,      1,
      "2",      "WEEK 1",          2, "DEF",       10,      2,
      "2",      "WEEK 2",          3, "DEF",       12,      3
    ) %>%
    mutate(STUDYID = "XYZ", .before = USUBJID)

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
      by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
      order = exprs(AVISITN, LBSEQ),
      mode = "last",
      set_values_to = exprs(
        AVISIT = "PBL LAST",
        AVISITN = 99,
        DTYPE = "LOV"
      )
    )
    #> # A tibble: 16 × 8
    #>    STUDYID USUBJID AVISIT   AVISITN PARAMCD  AVAL LBSEQ DTYPE
    #>    <chr>   <chr>   <chr>      <dbl> <chr>   <dbl> <dbl> <chr>
    #>  1 XYZ     1       BASELINE       1 ABC       120     1 <NA>
    #>  2 XYZ     1       WEEK 1         2 ABC       113     2 <NA>
    #>  3 XYZ     1       WEEK 1         2 ABC       117     3 <NA>
    #>  4 XYZ     2       BASELINE       1 ABC       101     1 <NA>
    #>  5 XYZ     2       WEEK 1         2 ABC       101     2 <NA>
    #>  6 XYZ     2       WEEK 2         3 ABC        95     3 <NA>
    #>  7 XYZ     1       BASELINE       1 DEF        17     1 <NA>
    #>  8 XYZ     1       WEEK 1         2 DEF        NA     2 <NA>
    #>  9 XYZ     1       WEEK 1         2 DEF        13     3 <NA>
    #> 10 XYZ     2       BASELINE       1 DEF         9     1 <NA>
    #> 11 XYZ     2       WEEK 1         2 DEF        10     2 <NA>
    #> 12 XYZ     2       WEEK 2         3 DEF        12     3 <NA>
    #> 13 XYZ     1       PBL LAST      99 ABC       117     3 LOV
    #> 14 XYZ     1       PBL LAST      99 DEF        13     3 LOV
    #> 15 XYZ     2       PBL LAST      99 ABC        95     3 LOV
    #> 16 XYZ     2       PBL LAST      99 DEF        12     3 LOV  

### Derive `DTYPE = "MINIMUM"`

For each subject and parameter the record with the minimum analysis
value should be selected and added as a new record to the input dataset.
If there are multiple records meeting the minimum, the first record with
respect to `AVISIT` and `LBSEQ` should be selected. For the new records
set `AVISIT = "PBL MIN"`, `AVISITN = 97`, and `DTYPE = "MINIMUM"`.

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
      by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
      order = exprs(AVAL, AVISITN, LBSEQ),
      mode = "first",
      set_values_to = exprs(
        AVISIT = "PBL MIN",
        AVISITN = 97,
        DTYPE = "MINIMUM"
      )
    )
    #> # A tibble: 16 × 8
    #>    STUDYID USUBJID AVISIT   AVISITN PARAMCD  AVAL LBSEQ DTYPE
    #>    <chr>   <chr>   <chr>      <dbl> <chr>   <dbl> <dbl> <chr>
    #>  1 XYZ     1       BASELINE       1 ABC       120     1 <NA>
    #>  2 XYZ     1       WEEK 1         2 ABC       113     2 <NA>
    #>  3 XYZ     1       WEEK 1         2 ABC       117     3 <NA>
    #>  4 XYZ     2       BASELINE       1 ABC       101     1 <NA>
    #>  5 XYZ     2       WEEK 1         2 ABC       101     2 <NA>
    #>  6 XYZ     2       WEEK 2         3 ABC        95     3 <NA>
    #>  7 XYZ     1       BASELINE       1 DEF        17     1 <NA>
    #>  8 XYZ     1       WEEK 1         2 DEF        NA     2 <NA>
    #>  9 XYZ     1       WEEK 1         2 DEF        13     3 <NA>
    #> 10 XYZ     2       BASELINE       1 DEF         9     1 <NA>
    #> 11 XYZ     2       WEEK 1         2 DEF        10     2 <NA>
    #> 12 XYZ     2       WEEK 2         3 DEF        12     3 <NA>
    #> 13 XYZ     1       PBL MIN       97 ABC       113     2 MINIMUM
    #> 14 XYZ     1       PBL MIN       97 DEF        13     3 MINIMUM
    #> 15 XYZ     2       PBL MIN       97 ABC        95     3 MINIMUM
    #> 16 XYZ     2       PBL MIN       97 DEF        10     2 MINIMUM

### Derive `DTYPE = "MAXIMUM"`

For each subject and parameter the record with the maximum analysis
value should be selected and added as a new record to the input dataset.
If there are multiple records meeting the maximum, the first record with
respect to `AVISIT` and `LBSEQ` should be selected. For the new records
set `AVISIT = "PBL MAX"`, `AVISITN = 98`, and `DTYPE = "MAXIMUM"`.

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
      by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
      order = exprs(desc(AVAL), AVISITN, LBSEQ),
      mode = "first",
      set_values_to = exprs(
        AVISIT = "PBL MAX",
        AVISITN = 99,
        DTYPE = "MAXIMUM"
      )
    )
    #> # A tibble: 16 × 8
    #>    STUDYID USUBJID AVISIT   AVISITN PARAMCD  AVAL LBSEQ DTYPE
    #>    <chr>   <chr>   <chr>      <dbl> <chr>   <dbl> <dbl> <chr>
    #>  1 XYZ     1       BASELINE       1 ABC       120     1 <NA>
    #>  2 XYZ     1       WEEK 1         2 ABC       113     2 <NA>
    #>  3 XYZ     1       WEEK 1         2 ABC       117     3 <NA>
    #>  4 XYZ     2       BASELINE       1 ABC       101     1 <NA>
    #>  5 XYZ     2       WEEK 1         2 ABC       101     2 <NA>
    #>  6 XYZ     2       WEEK 2         3 ABC        95     3 <NA>
    #>  7 XYZ     1       BASELINE       1 DEF        17     1 <NA>
    #>  8 XYZ     1       WEEK 1         2 DEF        NA     2 <NA>
    #>  9 XYZ     1       WEEK 1         2 DEF        13     3 <NA>
    #> 10 XYZ     2       BASELINE       1 DEF         9     1 <NA>
    #> 11 XYZ     2       WEEK 1         2 DEF        10     2 <NA>
    #> 12 XYZ     2       WEEK 2         3 DEF        12     3 <NA>
    #> 13 XYZ     1       PBL MAX       99 ABC       117     3 MAXIMUM
    #> 14 XYZ     1       PBL MAX       99 DEF        13     3 MAXIMUM
    #> 15 XYZ     2       PBL MAX       99 ABC       101     2 MAXIMUM
    #> 16 XYZ     2       PBL MAX       99 DEF        12     3 MAXIMUM

### Derive `DTYPE = "WOC"` or `DTYPE = "BOC"`

For each subject and parameter the record with the worst analysis value
should be selected and added as a new record to the input dataset. The
worst value is either the minimum or maximum value depending on the
parameter. If there are multiple records meeting the worst value, the
first record with respect to `AVISIT` and `LBSEQ` should be selected.
For the new records set `AVISIT = "PBL WORST"`, `AVISITN = 96`, and
`DTYPE = "WOC"`.

Here the maximum is considered worst for `PARAMCD = "ABC"` and the
minimum for `PARAMCD = "DEF"`.

    derive_extreme_records(
      adlb,
      dataset_add = adlb,
      filter_add = !is.na(AVAL) & AVISIT != "BASELINE",
      by_vars = exprs(!!!get_admiral_option("subject_keys"), PARAMCD),
      order = exprs(
        if_else(PARAMCD == "ABC", desc(AVAL), AVAL),
        AVISITN, LBSEQ
      ),
      mode = "first",
      set_values_to = exprs(
        AVISIT = "PBL WORST",
        AVISITN = 96,
        DTYPE = "WOC"
      )
    )
    #> # A tibble: 16 × 8
    #>    STUDYID USUBJID AVISIT    AVISITN PARAMCD  AVAL LBSEQ DTYPE
    #>    <chr>   <chr>   <chr>       <dbl> <chr>   <dbl> <dbl> <chr>
    #>  1 XYZ     1       BASELINE        1 ABC       120     1 <NA>
    #>  2 XYZ     1       WEEK 1          2 ABC       113     2 <NA>
    #>  3 XYZ     1       WEEK 1          2 ABC       117     3 <NA>
    #>  4 XYZ     2       BASELINE        1 ABC       101     1 <NA>
    #>  5 XYZ     2       WEEK 1          2 ABC       101     2 <NA>
    #>  6 XYZ     2       WEEK 2          3 ABC        95     3 <NA>
    #>  7 XYZ     1       BASELINE        1 DEF        17     1 <NA>
    #>  8 XYZ     1       WEEK 1          2 DEF        NA     2 <NA>
    #>  9 XYZ     1       WEEK 1          2 DEF        13     3 <NA>
    #> 10 XYZ     2       BASELINE        1 DEF         9     1 <NA>
    #> 11 XYZ     2       WEEK 1          2 DEF        10     2 <NA>
    #> 12 XYZ     2       WEEK 2          3 DEF        12     3 <NA>
    #> 13 XYZ     1       PBL WORST      96 ABC       117     3 WOC
    #> 14 XYZ     1       PBL WORST      96 DEF        13     3 WOC
    #> 15 XYZ     2       PBL WORST      96 ABC       101     2 WOC
    #> 16 XYZ     2       PBL WORST      96 DEF        10     2 WOC  

### Derive a parameter for the first disease progression (PD)

For each subject in the `ADSL` dataset a new parameter should be added
to the input dataset which indicates whether disease progression (PD)
occurred (set `AVALC = "Y"`, `AVAL = 1`) or not (set `AVALC = "N"`,
`AVAL = 0`). For the new parameter set `PARAMCD = "PD"` and
`PARAM = "Disease Progression"`.

    adsl <- tribble(
      ~USUBJID, ~DTHDT,
      "1",      ymd("2022-05-13"),
      "2",      ymd(""),
      "3",      ymd("")
    ) %>%
      mutate(STUDYID = "XX1234")

    adrs <- tribble(
      ~USUBJID, ~RSDTC,       ~AVALC, ~AVAL,
      "1",      "2020-01-02", "PR",       2,
      "1",      "2020-02-01", "CR",       1,
      "1",      "2020-03-01", "CR",       1,
      "2",      "2021-06-15", "SD",       3,
      "2",      "2021-07-16", "PD",       4,
      "2",      "2021-09-14", "PD",       4
    ) %>%
      mutate(
        STUDYID = "XX1234", .before = USUBJID
      ) %>%
      mutate(
        ADT = ymd(RSDTC),
        PARAMCD = "OVR",
        PARAM = "Overall Response",
        .after = RSDTC
      )

    derive_extreme_records(
      adrs,
      dataset_ref = adsl,
      dataset_add = adrs,
      by_vars = get_admiral_option("subject_keys"),
      filter_add = PARAMCD == "OVR" & AVALC == "PD",
      order = exprs(ADT),
      exist_flag = AVALC,
      true_value = "Y",
      false_value = "N",
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "PD",
        PARAM = "Disease Progression",
        AVAL = yn_to_numeric(AVALC),
      )
    )
    #> # A tibble: 9 × 8
    #>   STUDYID USUBJID RSDTC      ADT        PARAMCD PARAM               AVALC  AVAL
    #>   <chr>   <chr>   <chr>      <date>     <chr>   <chr>               <chr> <dbl>
    #> 1 XX1234  1       2020-01-02 2020-01-02 OVR     Overall Response    PR        2
    #> 2 XX1234  1       2020-02-01 2020-02-01 OVR     Overall Response    CR        1
    #> 3 XX1234  1       2020-03-01 2020-03-01 OVR     Overall Response    CR        1
    #> 4 XX1234  2       2021-06-15 2021-06-15 OVR     Overall Response    SD        3
    #> 5 XX1234  2       2021-07-16 2021-07-16 OVR     Overall Response    PD        4
    #> 6 XX1234  2       2021-09-14 2021-09-14 OVR     Overall Response    PD        4
    #> 7 XX1234  2       2021-07-16 2021-07-16 PD      Disease Progression Y         1
    #> 8 XX1234  1       <NA>       NA         PD      Disease Progression N         0
    #> 9 XX1234  3       <NA>       NA         PD      Disease Progression N         0

### Derive parameter indicating death

For each subject in the `ADSL` dataset a new parameter should be created
which indicates whether the subject died (set `AVALC = "Y"`, `AVAL = 1`)
or not (set `AVALC = "N"`, `AVAL = 0`). For the new parameter set
`PARAMCD = "DEATH"`, `PARAM = "Death"`, and `ADT` to the date of death
(`DTHDT`).

    derive_extreme_records(
      dataset_ref = adsl,
      dataset_add = adsl,
      by_vars = exprs(STUDYID, USUBJID),
      filter_add = !is.na(DTHDT),
      exist_flag = AVALC,
      true_value = "Y",
      false_value = "N",
      mode = "first",
      keep_source_vars = exprs(AVALC),
      set_values_to = exprs(
        PARAMCD = "DEATH",
        PARAM = "Death",
        ADT = DTHDT
      )
    )
    #> # A tibble: 3 × 6
    #>   STUDYID USUBJID PARAMCD PARAM ADT        AVALC
    #>   <chr>   <chr>   <chr>   <chr> <date>     <chr>
    #> 1 XX1234  1       DEATH   Death 2022-05-13 Y
    #> 2 XX1234  2       DEATH   Death NA         N
    #> 3 XX1234  3       DEATH   Death NA         N    

The `keep_source_vars` argument is specified to avoid that all `ADSL`
variables (like `DTHDT`) are copied to the parameter.
