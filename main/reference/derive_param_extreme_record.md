# Adds a Parameter Based on First or Last Record from Multiple Sources

**\[deprecated\]** The `derive_param_extreme_record()` function has been
deprecated in favor of
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md).

Generates parameter based on the first or last observation from multiple
source datasets, based on user-defined filter, order and by group
criteria. All variables of the selected observation are kept.

## Usage

``` r
derive_param_extreme_record(
  dataset = NULL,
  sources,
  source_datasets,
  by_vars = NULL,
  order,
  mode,
  set_values_to
)
```

## Arguments

- dataset:

  Input dataset

  Default value

  :   `NULL`

- sources:

  Sources

  A list of
  [`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md)
  objects is expected.

  Default value

  :   none

- source_datasets:

  Source datasets

  A named list of datasets is expected. The `dataset_name` field of
  [`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md)
  refers to the dataset provided in the list. The variables specified by
  the `order` and the `by_vars` arguments are expected after applying
  `new_vars`.

  Default value

  :   none

- by_vars:

  Grouping variables

  If the argument is specified, for each by group the observations are
  selected separately.

  Default value

  :   `NULL`

- order:

  Sort order

  If the argument is set to a non-null value, for each by group the
  first or last observation from the source datasets is selected with
  respect to the specified order. Variables created via `new_vars` e.g.,
  imputed date variables, can be specified as well (see examples below).

  Please note that `NA` is considered as the last value. I.e., if a
  order variable is `NA` and `mode = "last"`, this observation is chosen
  while for `mode = "first"` the observation is chosen only if there are
  no observations where the variable is not `NA`.

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))`

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, for each by group the first observation
  with respect to `order` is included in the output dataset. If `"last"`
  is specified, the last observation is included in the output dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations.

  A list of variable name-value pairs is expected.

  - LHS refers to a variable.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value or `NA`, e.g.,
    `exprs(PARAMCD = "PD", PARAM = "First Progressive Disease")`.

  Default value

  :   none

## Value

The input dataset with the first or last observation of each by group
added as new observations.

## Details

The following steps are performed to create the output dataset:

1.  For each source dataset the observations as specified by the
    `filter` element are selected.

2.  Variables specified by `new_vars` are created for each source
    dataset.

3.  The first or last observation (with respect to the `order` variable)
    for each by group (specified by `by_vars`) from multiple sources is
    selected and added to the input dataset.

## See also

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)

## Examples

``` r
aevent_samp <- tibble::tribble(
  ~USUBJID, ~PARAMCD,                       ~PARAM,     ~RSSTDTC,
  "1",          "PD",  "First Progressive Disease", "2022-04-01",
  "2",          "PD",  "First Progressive Disease", "2021-04-01",
  "3",          "PD",  "First Progressive Disease", "2023-04-01"
)

cm <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~CMDECOD,     ~CMSTDTC,
  "1001",        "1",    "ACT", "2021-12-25"
)

pr <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
  "1001",        "1",    "ACS", "2021-12-27",
  "1001",        "2",    "ACS", "2020-12-25",
  "1001",        "3",    "ACS", "2022-12-25",
)
derive_param_extreme_record(
  dataset = aevent_samp,
  sources = list(
    records_source(
      dataset_name = "cm",
      filter = CMDECOD == "ACT",
      new_vars = exprs(
        ADT = convert_dtc_to_dt(CMSTDTC),
        AVALC = CMDECOD
      )
    ),
    records_source(
      dataset_name = "pr",
      filter = PRDECOD == "ACS",
      new_vars = exprs(
        ADT = convert_dtc_to_dt(PRSTDTC),
        AVALC = PRDECOD
      )
    )
  ),
  source_datasets = list(cm = cm, pr = pr),
  by_vars = exprs(USUBJID),
  order = exprs(ADT),
  mode = "first",
  set_values_to = exprs(
    PARAMCD = "FIRSTACT",
    PARAM = "First Anti-Cancer Therapy"
  )
)
#> Warning: `derive_param_extreme_record()` was deprecated in admiral 1.2.0.
#> ℹ Please use `derive_extreme_event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
#> # A tibble: 6 × 6
#>   USUBJID PARAMCD  PARAM                     RSSTDTC    ADT        AVALC
#>   <chr>   <chr>    <chr>                     <chr>      <date>     <chr>
#> 1 1       PD       First Progressive Disease 2022-04-01 NA         NA   
#> 2 2       PD       First Progressive Disease 2021-04-01 NA         NA   
#> 3 3       PD       First Progressive Disease 2023-04-01 NA         NA   
#> 4 1       FIRSTACT First Anti-Cancer Therapy NA         2021-12-25 ACT  
#> 5 2       FIRSTACT First Anti-Cancer Therapy NA         2020-12-25 ACS  
#> 6 3       FIRSTACT First Anti-Cancer Therapy NA         2022-12-25 ACS  
```
