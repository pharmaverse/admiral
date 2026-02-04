# Consolidate Multiple Meta Datasets Into a Single One

The purpose of the function is to consolidate multiple meta datasets
into a single one. For example, from global and project specific
parameter mappings a single lookup table can be created.

## Usage

``` r
consolidate_metadata(
  datasets,
  key_vars,
  source_var = SOURCE,
  check_vars = "warning",
  check_type = "error"
)
```

## Arguments

- datasets:

  List of datasets to consolidate

  Permitted values

  :   A named list of datasets

  Default value

  :   none

- key_vars:

  Key variables

  The specified variables must be a unique of all input datasets.

  Permitted values

  :   A list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)

  Default value

  :   none

- source_var:

  Source variable

  The specified variable is added to the output dataset. It is set the
  name of the dataset the observation is originating from.

  Permitted values

  :   A symbol

  Default value

  :   `SOURCE`

- check_vars:

  Check variables?

  If `"message"`, `"warning"`, or `"error"` is specified, a message is
  issued if the variable names differ across the input datasets
  (`datasets`).

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, a message is issued if the
  key variables (`key_vars`) are not a unique key in all of the input
  datasets (`datasets`).

  Permitted values

  :   `"none"`, `"warning"`, `"error"`

  Default value

  :   `"error"`

## Value

A dataset which contains one row for each by group occurring in any of
the input datasets.

## Details

All observations of the input datasets are put together into a single
dataset. If a by group (defined by `key_vars`) exists in more than one
of the input datasets, the observation from the last dataset is
selected.

## See also

Creating auxiliary datasets:
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_period_dataset.md),
[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md),
[`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)

## Examples

``` r
library(tibble)
glob_ranges <- tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "PULSE",      60,    100,
  "SYSBP",      90,    130,
  "DIABP",      60,     80
)
proj_ranges <- tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "SYSBP",     100,    140,
  "DIABP",      70,     90
)
stud_ranges <- tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "BMI",        18,     25
)

consolidate_metadata(
  datasets = list(
    global = glob_ranges,
    project = proj_ranges,
    study = stud_ranges
  ),
  key_vars = exprs(PARAMCD)
)
#> # A tibble: 4 × 4
#>   SOURCE  PARAMCD ANRLO ANRHI
#>   <chr>   <chr>   <dbl> <dbl>
#> 1 study   BMI        18    25
#> 2 project DIABP      70    90
#> 3 global  PULSE      60   100
#> 4 project SYSBP     100   140
```
