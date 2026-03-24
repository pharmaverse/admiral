# Create Summary Records

**\[deprecated\]** The `get_summary_records()` has been deprecated in
favor of
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
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
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
function.

## See also

[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)

## Examples

``` r
library(tibble)

adeg <- tribble(
  ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
  "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
  "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
  "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
  "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:45", 384,   "Placebo",
  "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
  "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
  "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:45", 385,   "Placebo",
  "XYZ-1001", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
  "XYZ-1001", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
  "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
  "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
  "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
  "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:50", 401,   "Active 20mg",
  "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
  "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
  "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:50", 412,   "Active 20mg",
  "XYZ-1002", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
  "XYZ-1002", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
)

# Summarize the average of the triplicate ECG interval values (AVAL)
get_summary_records(
  adeg,
  by_vars = exprs(USUBJID, PARAM, AVISIT),
  set_values_to = exprs(
    AVAL = mean(AVAL, na.rm = TRUE),
    DTYPE = "AVERAGE"
  )
)
#> Warning: `get_summary_records()` was deprecated in admiral 1.2.0.
#> ℹ Please use `derive_summary_records()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
#> # A tibble: 6 × 5
#>   USUBJID  PARAM            AVISIT    AVAL DTYPE  
#>   <chr>    <chr>            <chr>    <dbl> <chr>  
#> 1 XYZ-1001 QTcF Int. (msec) Baseline  393. AVERAGE
#> 2 XYZ-1001 QTcF Int. (msec) Visit 2   388. AVERAGE
#> 3 XYZ-1001 QTcF Int. (msec) Visit 3   394. AVERAGE
#> 4 XYZ-1002 QTcF Int. (msec) Baseline  400. AVERAGE
#> 5 XYZ-1002 QTcF Int. (msec) Visit 2   403. AVERAGE
#> 6 XYZ-1002 QTcF Int. (msec) Visit 3   409. AVERAGE

# Derive more than one summary variable
get_summary_records(
  adeg,
  by_vars = exprs(USUBJID, PARAM, AVISIT),
  set_values_to = exprs(
    AVAL = mean(AVAL),
    ASTDTM = min(convert_dtc_to_dtm(EGDTC)),
    AENDTM = max(convert_dtc_to_dtm(EGDTC)),
    DTYPE = "AVERAGE"
  )
)
#> # A tibble: 6 × 7
#>   USUBJID  PARAM      AVISIT  AVAL ASTDTM              AENDTM              DTYPE
#>   <chr>    <chr>      <chr>  <dbl> <dttm>              <dttm>              <chr>
#> 1 XYZ-1001 QTcF Int.… Basel…  393. 2016-02-24 07:50:00 2016-02-24 07:56:00 AVER…
#> 2 XYZ-1001 QTcF Int.… Visit…  388. 2016-03-08 09:45:00 2016-03-08 09:51:00 AVER…
#> 3 XYZ-1001 QTcF Int.… Visit…  394. 2016-03-22 10:45:00 2016-03-22 10:51:00 AVER…
#> 4 XYZ-1002 QTcF Int.… Basel…  400. 2016-02-22 07:58:00 2016-02-22 08:01:00 AVER…
#> 5 XYZ-1002 QTcF Int.… Visit…  403. 2016-03-06 09:50:00 2016-03-06 09:56:00 AVER…
#> 6 XYZ-1002 QTcF Int.… Visit…  409. 2016-03-24 10:50:00 2016-03-24 10:56:00 AVER…

# Sample ADEG dataset with triplicate record for only AVISIT = 'Baseline'
adeg <- tribble(
  ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
  "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
  "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
  "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
  "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
  "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
  "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
  "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
  "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
  "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
  "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
  "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
  "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
  "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
  "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
)

# Compute the average of AVAL only if there are more than 2 records within the
# by group
get_summary_records(
  adeg,
  by_vars = exprs(USUBJID, PARAM, AVISIT),
  filter = n() > 2,
  set_values_to = exprs(
    AVAL = mean(AVAL, na.rm = TRUE),
    DTYPE = "AVERAGE"
  )
)
#> # A tibble: 2 × 5
#>   USUBJID  PARAM            AVISIT    AVAL DTYPE  
#>   <chr>    <chr>            <chr>    <dbl> <chr>  
#> 1 XYZ-1001 QTcF Int. (msec) Baseline  393. AVERAGE
#> 2 XYZ-1002 QTcF Int. (msec) Baseline  400. AVERAGE
```
