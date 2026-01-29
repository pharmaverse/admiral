# Filter the First or Last Observation for Each By Group

Filters the first or last observation for each by group.

## Usage

``` r
filter_extreme(dataset, by_vars = NULL, order, mode, check_type = "warning")
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- order:

  Sort order

  Within each by group the observations are ordered by the specified
  order.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, the first observation of each by group is
  included in the output dataset. If `"last"` is specified, the last
  observation of each by group is included in the output dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the input dataset are not unique with
  respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

A dataset containing the first or last observation of each by group

## Details

For each group (with respect to the variables specified for the
`by_vars` parameter) the first or last observation (with respect to the
order specified for the `order` parameter and the mode specified for the
`mode` parameter) is included in the output dataset.

## See also

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md),
[`filter_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_exist.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_not_exist.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/main/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)

ex <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~EXSEQ, ~EXDOSE,    ~EXTRT,
  "PILOT01",    "EX", "01-1442",      1,      54,    "XANO",
  "PILOT01",    "EX", "01-1442",      2,      54,    "XANO",
  "PILOT01",    "EX", "01-1442",      3,      54,    "XANO",
  "PILOT01",    "EX", "01-1444",      1,      54,    "XANO",
  "PILOT01",    "EX", "01-1444",      2,      81,    "XANO",
  "PILOT01",    "EX", "05-1382",      1,      54,    "XANO",
  "PILOT01",    "EX", "08-1213",      1,      54,    "XANO",
  "PILOT01",    "EX", "10-1053",      1,      54,    "XANO",
  "PILOT01",    "EX", "10-1053",      2,      54,    "XANO",
  "PILOT01",    "EX", "10-1183",      1,       0, "PLACEBO",
  "PILOT01",    "EX", "10-1183",      2,       0, "PLACEBO",
  "PILOT01",    "EX", "10-1183",      3,       0, "PLACEBO",
  "PILOT01",    "EX", "11-1036",      1,       0, "PLACEBO",
  "PILOT01",    "EX", "11-1036",      2,       0, "PLACEBO",
  "PILOT01",    "EX", "11-1036",      3,       0, "PLACEBO",
  "PILOT01",    "EX", "14-1425",      1,      54,    "XANO",
  "PILOT01",    "EX", "15-1319",      1,      54,    "XANO",
  "PILOT01",    "EX", "15-1319",      2,      81,    "XANO",
  "PILOT01",    "EX", "16-1151",      1,      54,    "XANO",
  "PILOT01",    "EX", "16-1151",      2,      54,    "XANO"
)


# Select first dose for each patient
ex %>%
  filter_extreme(
    by_vars = exprs(USUBJID),
    order = exprs(EXSEQ),
    mode = "first"
  ) %>%
  select(USUBJID, EXSEQ)
#> # A tibble: 10 × 2
#>    USUBJID EXSEQ
#>    <chr>   <dbl>
#>  1 01-1442     1
#>  2 01-1444     1
#>  3 05-1382     1
#>  4 08-1213     1
#>  5 10-1053     1
#>  6 10-1183     1
#>  7 11-1036     1
#>  8 14-1425     1
#>  9 15-1319     1
#> 10 16-1151     1

# Select highest dose for each patient on the active drug
ex %>%
  filter(EXTRT != "PLACEBO") %>%
  filter_extreme(
    by_vars = exprs(USUBJID),
    order = exprs(EXDOSE),
    mode = "last",
    check_type = "none"
  ) %>%
  select(USUBJID, EXTRT, EXDOSE)
#> # A tibble: 8 × 3
#>   USUBJID EXTRT EXDOSE
#>   <chr>   <chr>  <dbl>
#> 1 01-1442 XANO      54
#> 2 01-1444 XANO      81
#> 3 05-1382 XANO      54
#> 4 08-1213 XANO      54
#> 5 10-1053 XANO      54
#> 6 14-1425 XANO      54
#> 7 15-1319 XANO      81
#> 8 16-1151 XANO      54
```
