# Derive LOCF (Last Observation Carried Forward) Records

Adds LOCF records as new observations for each 'by group' when the
dataset does not contain observations for missed visits/time points and
when analysis value is missing.

## Usage

``` r
derive_locf_records(
  dataset,
  dataset_ref,
  by_vars,
  id_vars_ref = NULL,
  analysis_var = AVAL,
  imputation = "add",
  order,
  keep_vars = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars`, `analysis_var`, `order`, and
  `keep_vars` arguments are expected to be in the dataset.

  Default value

  :   none

- dataset_ref:

  Expected observations dataset

  Data frame with all the combinations of `PARAMCD`, `PARAM`, `AVISIT`,
  `AVISITN`, ... which are expected in the dataset is expected.

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` those observations from
  `dataset_ref` are added to the output dataset which do not have a
  corresponding observation in the input dataset or for which
  `analysis_var` is `NA` for the corresponding observation in the input
  dataset.

  Default value

  :   none

- id_vars_ref:

  Grouping variables in expected observations dataset

  The variables to group by in `dataset_ref` when determining which
  observations should be added to the input dataset.

  Default value

  :   All the variables in `dataset_ref`

- analysis_var:

  Analysis variable.

  Permitted values

  :   a variable

  Default value

  :   `AVAL`

- imputation:

  Select the mode of imputation:

  `add`: Keep all original records and add imputed records for missing
  timepoints and missing `analysis_var` values from `dataset_ref`.

  `update`: Update records with missing `analysis_var` and add imputed
  records for missing timepoints from `dataset_ref`.

  `update_add`: Keep all original records, update records with missing
  `analysis_var` and add imputed records for missing timepoints from
  `dataset_ref`.

  Permitted values

  :   One of these 3 values: `"add"`, `"update"`, `"update_add"`

  Default value

  :   `"add"`

- order:

  Sort order

  The dataset is sorted by `order` before carrying the last observation
  forward (e.g. `AVAL`) within each `by_vars`.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/test_cicd/articles/generic.md).

  Default value

  :   none

- keep_vars:

  Variables that need carrying the last observation forward

  Keep variables that need carrying the last observation forward other
  than `analysis_var` (e.g., `PARAMN`, `VISITNUM`). If by default
  `NULL`, only variables specified in `by_vars` and `analysis_var` will
  be populated in the newly created records.

  Default value

  :   `NULL`

## Value

The input dataset with the new "LOCF" observations added for each
`by_vars`, based on the value passed to the `imputation` argument.

## Details

For each group (with respect to the variables specified for the by_vars
parameter) those observations from `dataset_ref` are added to the output
dataset

- which do not have a corresponding observation in the input dataset or

- for which `analysis_var` is `NA` for the corresponding observation in
  the input dataset.

  For the new observations, `analysis_var` is set to the non-missing
  `analysis_var` of the previous observation in the input dataset (when
  sorted by `order`) and `DTYPE` is set to "LOCF".

  The `imputation` argument decides whether to update the existing
  observation when `analysis_var` is `NA` (`"update"` and
  `"update_add"`), or to add a new observation from `dataset_ref`
  instead (`"add"`).

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_records.md),
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

## Author

G Gayatri

## Examples

### Add records for missing analysis variable using reference dataset

Imputed records should be added for missing timepoints and for missing
`analysis_var` (from `dataset_ref`), while retaining all original
records.

- The reference dataset for the imputed records is specified by the
  `dataset_add` argument. It should contain all expected combinations of
  variables. In this case, `advs_expected_obsv` is created by
  `crossing()` datasets `paramcd` and `avisit`, which includes all
  combinations of PARAMCD, AVISITN, and AVISIT.

- The groups for which new records are added are specified by the
  `by_vars` argument. Here, one record should be added for each
  *subject* and *parameter*. Therefore,
  `by_vars = exprs(STUDYID, USUBJID, PARAMCD)` is specified.

- The imputation method is specified using the `imputation` argument. In
  this case, records with missing analysis values *add* records from
  `dataset_ref` after the data are sorted by the variables in `by_vars`
  and by visit (`AVISITN` and `AVISIT`), as specified in the `order`
  argument.

- Variables other than `analysis_var` and `by_vars` that require LOCF
  (Last-Observation- Carried-Forward handling (in this case, `PARAMN`)
  are specified in the `keep_vars` argument.

    library(dplyr)
    library(tibble)
    library(tidyr)

    advs <- tribble(
      ~STUDYID,  ~USUBJID,      ~VSSEQ, ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT,
      "CDISC01", "01-701-1015",      1, "PULSE",        1,    65,        0, "BASELINE",
      "CDISC01", "01-701-1015",      2, "DIABP",        2,    79,        0, "BASELINE",
      "CDISC01", "01-701-1015",      3, "DIABP",        2,    80,        2, "WEEK 2",
      "CDISC01", "01-701-1015",      4, "DIABP",        2,    NA,        4, "WEEK 4",
      "CDISC01", "01-701-1015",      5, "DIABP",        2,    NA,        6, "WEEK 6",
      "CDISC01", "01-701-1015",      6, "SYSBP",        3,   130,        0, "BASELINE",
      "CDISC01", "01-701-1015",      7, "SYSBP",        3,   132,        2, "WEEK 2"
    )

    paramcd <- tribble(
    ~PARAMCD,
    "PULSE",
    "DIABP",
    "SYSBP"
    )

    avisit <- tribble(
    ~AVISITN, ~AVISIT,
    0, "BASELINE",
    2, "WEEK 2",
    4, "WEEK 4",
    6, "WEEK 6"
    )

    advs_expected_obsv <- paramcd %>%
    crossing(avisit)

    derive_locf_records(
      dataset = advs,
      dataset_ref = advs_expected_obsv,
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      imputation = "add",
      order = exprs(AVISITN, AVISIT),
      keep_vars = exprs(PARAMN)
    ) |>
      arrange(USUBJID, PARAMCD, AVISIT)
    #> # A tibble: 14 × 9
    #>    STUDYID USUBJID     VSSEQ PARAMCD PARAMN  AVAL AVISITN AVISIT   DTYPE
    #>    <chr>   <chr>       <dbl> <chr>    <dbl> <dbl>   <dbl> <chr>    <chr>
    #>  1 CDISC01 01-701-1015     2 DIABP        2    79       0 BASELINE <NA>
    #>  2 CDISC01 01-701-1015     3 DIABP        2    80       2 WEEK 2   <NA>
    #>  3 CDISC01 01-701-1015    NA DIABP        2    80       4 WEEK 4   LOCF
    #>  4 CDISC01 01-701-1015     4 DIABP        2    NA       4 WEEK 4   <NA>
    #>  5 CDISC01 01-701-1015    NA DIABP        2    80       6 WEEK 6   LOCF
    #>  6 CDISC01 01-701-1015     5 DIABP        2    NA       6 WEEK 6   <NA>
    #>  7 CDISC01 01-701-1015     1 PULSE        1    65       0 BASELINE <NA>
    #>  8 CDISC01 01-701-1015    NA PULSE        1    65       2 WEEK 2   LOCF
    #>  9 CDISC01 01-701-1015    NA PULSE        1    65       4 WEEK 4   LOCF
    #> 10 CDISC01 01-701-1015    NA PULSE        1    65       6 WEEK 6   LOCF
    #> 11 CDISC01 01-701-1015     6 SYSBP        3   130       0 BASELINE <NA>
    #> 12 CDISC01 01-701-1015     7 SYSBP        3   132       2 WEEK 2   <NA>
    #> 13 CDISC01 01-701-1015    NA SYSBP        3   132       4 WEEK 4   LOCF
    #> 14 CDISC01 01-701-1015    NA SYSBP        3   132       6 WEEK 6   LOCF 

### Update records for missing analysis variable

When the `imputation` mode is set to *update*, missing `analysis_var`
values are updated using values from the last record after the dataset
is sorted by `by_vars` and `order`. Imputed records are added for
missing timepoints (from `dataset_ref`).

    derive_locf_records(
      dataset = advs,
      dataset_ref = advs_expected_obsv,
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      imputation = "update",
      order = exprs(AVISITN, AVISIT),
    ) |>
      arrange(USUBJID, PARAMCD, AVISIT)
    #> # A tibble: 12 × 9
    #>    STUDYID USUBJID     VSSEQ PARAMCD PARAMN  AVAL AVISITN AVISIT   DTYPE
    #>    <chr>   <chr>       <dbl> <chr>    <dbl> <dbl>   <dbl> <chr>    <chr>
    #>  1 CDISC01 01-701-1015     2 DIABP        2    79       0 BASELINE <NA>
    #>  2 CDISC01 01-701-1015     3 DIABP        2    80       2 WEEK 2   <NA>
    #>  3 CDISC01 01-701-1015     4 DIABP        2    80       4 WEEK 4   LOCF
    #>  4 CDISC01 01-701-1015     5 DIABP        2    80       6 WEEK 6   LOCF
    #>  5 CDISC01 01-701-1015     1 PULSE        1    65       0 BASELINE <NA>
    #>  6 CDISC01 01-701-1015    NA PULSE       NA    65       2 WEEK 2   LOCF
    #>  7 CDISC01 01-701-1015    NA PULSE       NA    65       4 WEEK 4   LOCF
    #>  8 CDISC01 01-701-1015    NA PULSE       NA    65       6 WEEK 6   LOCF
    #>  9 CDISC01 01-701-1015     6 SYSBP        3   130       0 BASELINE <NA>
    #> 10 CDISC01 01-701-1015     7 SYSBP        3   132       2 WEEK 2   <NA>
    #> 11 CDISC01 01-701-1015    NA SYSBP       NA   132       4 WEEK 4   LOCF
    #> 12 CDISC01 01-701-1015    NA SYSBP       NA   132       6 WEEK 6   LOCF 

### Update records for missing analysis variable while keeping the original records

When the `imputation` mode is set to *update_add*, the missing
`analysis_var` values are updated using values from the last record
after the dataset is sorted by `by_vars` and `order`. The updated values
are added as new records, while the original records with missing
`analysis_var` are retained. Imputed records are added for missing
timepoints (from `dataset_ref`).

    derive_locf_records(
      dataset = advs,
      dataset_ref = advs_expected_obsv,
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      imputation = "update_add",
      order = exprs(AVISITN, AVISIT),
    ) |>
      arrange(USUBJID, PARAMCD, AVISIT)
    #> # A tibble: 14 × 9
    #>    STUDYID USUBJID     VSSEQ PARAMCD PARAMN  AVAL AVISITN AVISIT   DTYPE
    #>    <chr>   <chr>       <dbl> <chr>    <dbl> <dbl>   <dbl> <chr>    <chr>
    #>  1 CDISC01 01-701-1015     2 DIABP        2    79       0 BASELINE <NA>
    #>  2 CDISC01 01-701-1015     3 DIABP        2    80       2 WEEK 2   <NA>
    #>  3 CDISC01 01-701-1015     4 DIABP        2    80       4 WEEK 4   LOCF
    #>  4 CDISC01 01-701-1015     4 DIABP        2    NA       4 WEEK 4   <NA>
    #>  5 CDISC01 01-701-1015     5 DIABP        2    80       6 WEEK 6   LOCF
    #>  6 CDISC01 01-701-1015     5 DIABP        2    NA       6 WEEK 6   <NA>
    #>  7 CDISC01 01-701-1015     1 PULSE        1    65       0 BASELINE <NA>
    #>  8 CDISC01 01-701-1015    NA PULSE       NA    65       2 WEEK 2   LOCF
    #>  9 CDISC01 01-701-1015    NA PULSE       NA    65       4 WEEK 4   LOCF
    #> 10 CDISC01 01-701-1015    NA PULSE       NA    65       6 WEEK 6   LOCF
    #> 11 CDISC01 01-701-1015     6 SYSBP        3   130       0 BASELINE <NA>
    #> 12 CDISC01 01-701-1015     7 SYSBP        3   132       2 WEEK 2   <NA>
    #> 13 CDISC01 01-701-1015    NA SYSBP       NA   132       4 WEEK 4   LOCF
    #> 14 CDISC01 01-701-1015    NA SYSBP       NA   132       6 WEEK 6   LOCF 
