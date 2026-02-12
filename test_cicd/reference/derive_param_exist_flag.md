# Add an Existence Flag Parameter

Add a new parameter indicating that a certain event exists in a dataset.
`AVALC` and `AVAL` indicate if an event occurred or not. For example,
the function can derive a parameter indicating if there is measurable
disease at baseline.

## Usage

``` r
derive_param_exist_flag(
  dataset = NULL,
  dataset_ref,
  dataset_add,
  condition,
  true_value = "Y",
  false_value = NA_character_,
  missing_value = NA_character_,
  filter_add = NULL,
  by_vars = get_admiral_option("subject_keys"),
  set_values_to
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD` is expected as well.

  Default value

  :   `NULL`

- dataset_ref:

  Reference dataset, e.g., ADSL

  The variables specified in `by_vars` are expected. For each group (as
  defined by `by_vars`) from the specified dataset (`dataset_ref`), the
  existence flag is calculated and added as a new observation to the
  input datasets (`dataset`).

  Default value

  :   none

- dataset_add:

  Additional dataset

  The variables specified by the `by_vars` parameter are expected.

  This dataset is used to check if an event occurred or not. Any
  observation in the dataset fulfilling the event condition
  (`condition`) is considered as an event.

  Default value

  :   none

- condition:

  Event condition

  The condition is evaluated at the additional dataset (`dataset_add`).

  For all groups where it evaluates as `TRUE` at least once `AVALC` is
  set to the true value (`true_value`) for the new observations.

  For all groups where it evaluates as `FALSE` or `NA` for all
  observations `AVALC` is set to the false value (`false_value`).

  For all groups not present in the additional dataset `AVALC` is set to
  the missing value (`missing_value`).

  Default value

  :   none

- true_value:

  True value

  For all groups with at least one observations in the additional
  dataset (`dataset_add`) fulfilling the event condition (`condition`),
  `AVALC` is set to the specified value (`true_value`).

  Permitted values

  :   A character scalar

  Default value

  :   `"Y"`

- false_value:

  False value

  For all groups with at least one observations in the additional
  dataset (`dataset_add`) but none of them is fulfilling the event
  condition (`condition`), `AVALC` is set to the specified value
  (`false_value`).

  Permitted values

  :   A character scalar

  Default value

  :   `NA_character_`

- missing_value:

  Values used for missing information

  For all groups without an observation in the additional dataset
  (`dataset_add`), `AVALC` is set to the specified value
  (`missing_value`).

  Permitted values

  :   A character scalar

  Default value

  :   `NA_character_`

- filter_add:

  Filter for additional data

  Only observations fulfilling the specified condition are taken into
  account for flagging. If the parameter is not specified, all
  observations are considered.

  Permitted values

  :   a condition

  Default value

  :   `NULL`

- by_vars:

  Grouping variables

  Default value

  :   `get_admiral_option("subject_keys")`

- set_values_to:

  Variables to set

  A named list returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  defining the variables to be set for the new parameter, e.g.
  `exprs(PARAMCD = "MDIS", PARAM = "Measurable Disease at Baseline")` is
  expected. The values must be symbols, character strings, numeric
  values, `NA`, or expressions.

  Default value

  :   none

## Value

The input dataset with a new parameter indicating if an event occurred
(`AVALC` and the variables specified by `by_vars` and `set_value_to` are
populated for the new parameter).

## Details

1.  The additional dataset (`dataset_add`) is restricted to the
    observations matching the `filter_add` condition.

2.  For each group in `dataset_ref` a new observation is created.

    - The `AVALC` variable is added and set to the true value
      (`true_value`) if for the group at least one observation exists in
      the (restricted) additional dataset where the condition evaluates
      to `TRUE`.

    - It is set to the false value (`false_value`) if for the group at
      least one observation exists and for all observations the
      condition evaluates to `FALSE` or `NA`.

    - Otherwise, it is set to the missing value (`missing_value`), i.e.,
      for those groups not in `dataset_add`.

3.  The variables specified by the `set_values_to` parameter are added
    to the new observations.

4.  The new observations are added to input dataset.

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_doseint.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

# Derive a new parameter for measurable disease at baseline
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
) %>%
  mutate(STUDYID = "XX1234")

tu <- tribble(
  ~USUBJID, ~VISIT,      ~TUSTRESC,
  "1",      "SCREENING", "TARGET",
  "1",      "WEEK 1",    "TARGET",
  "1",      "WEEK 5",    "TARGET",
  "1",      "WEEK 9",    "NON-TARGET",
  "2",      "SCREENING", "NON-TARGET",
  "2",      "SCREENING", "NON-TARGET"
) %>%
  mutate(
    STUDYID = "XX1234",
    TUTESTCD = "TUMIDENT"
  )

derive_param_exist_flag(
  dataset_ref = adsl,
  dataset_add = tu,
  filter_add = TUTESTCD == "TUMIDENT" & VISIT == "SCREENING",
  condition = TUSTRESC == "TARGET",
  false_value = "N",
  missing_value = "N",
  set_values_to = exprs(
    AVAL = yn_to_numeric(AVALC),
    PARAMCD = "MDIS",
    PARAM = "Measurable Disease at Baseline"
  )
)
#> # A tibble: 3 × 6
#>   STUDYID USUBJID AVALC  AVAL PARAMCD PARAM                         
#>   <chr>   <chr>   <chr> <dbl> <chr>   <chr>                         
#> 1 XX1234  1       Y         1 MDIS    Measurable Disease at Baseline
#> 2 XX1234  2       N         0 MDIS    Measurable Disease at Baseline
#> 3 XX1234  3       N         0 MDIS    Measurable Disease at Baseline
```
