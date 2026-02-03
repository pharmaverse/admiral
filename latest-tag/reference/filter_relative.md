# Filter the Observations Before or After a Condition is Fulfilled

Filters the observations before or after the observation where a
specified condition is fulfilled for each by group. For example, the
function could be called to select for each subject all observations
before the first disease progression.

## Usage

``` r
filter_relative(
  dataset,
  by_vars,
  order,
  condition,
  mode,
  selection,
  inclusive,
  keep_no_ref_groups = TRUE,
  check_type = "warning"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   none

- order:

  Sort order

  Within each by group the observations are ordered by the specified
  order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))`

  Default value

  :   none

- condition:

  Condition for Reference Observation

  The specified condition determines the reference observation. The
  output dataset contains all observations before or after (`selection`
  parameter) the reference observation.

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, for each by group the observations before
  or after (`selection` parameter) the observation where the condition
  (`condition` parameter) is fulfilled the *first* time is included in
  the output dataset. If `"last"` is specified, for each by group the
  observations before or after (`selection` parameter) the observation
  where the condition (`condition` parameter) is fulfilled the *last*
  time is included in the output dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- selection:

  Select observations before or after the reference observation?

  Permitted values

  :   `"before"`, `"after"`

  Default value

  :   none

- inclusive:

  Include the reference observation?

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   none

- keep_no_ref_groups:

  Should by groups without reference observation be kept?

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the input dataset are not unique with
  respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

A dataset containing for each by group the observations before or after
the observation where the condition was fulfilled the first or last time

## Details

For each by group ( `by_vars` parameter) the observations before or
after (`selection` parameter) the observations where the condition
(`condition` parameter) is fulfilled the first or last time (`order`
parameter and `mode` parameter) is included in the output dataset.

## See also

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/count_vals.md),
[`filter_exist()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_exist.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_extreme.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_joined.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_not_exist.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/min_cond.md)

## Examples

``` r
library(tibble)

response <- tribble(
  ~USUBJID, ~AVISITN, ~AVALC,
  "1",      1,        "PR",
  "1",      2,        "CR",
  "1",      3,        "CR",
  "1",      4,        "SD",
  "1",      5,        "NE",
  "2",      1,        "SD",
  "2",      2,        "PD",
  "2",      3,        "PD",
  "3",      1,        "SD",
  "4",      1,        "SD",
  "4",      2,        "PR",
  "4",      3,        "PD",
  "4",      4,        "SD",
  "4",      5,        "PR"
)

# Select observations up to first PD for each patient
response %>%
  filter_relative(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )
#> # A tibble: 11 × 3
#>    USUBJID AVISITN AVALC
#>    <chr>     <dbl> <chr>
#>  1 1             1 PR   
#>  2 1             2 CR   
#>  3 1             3 CR   
#>  4 1             4 SD   
#>  5 1             5 NE   
#>  6 2             1 SD   
#>  7 2             2 PD   
#>  8 3             1 SD   
#>  9 4             1 SD   
#> 10 4             2 PR   
#> 11 4             3 PD   

# Select observations after last CR, PR, or SD for each patient
response %>%
  filter_relative(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    condition = AVALC %in% c("CR", "PR", "SD"),
    mode = "last",
    selection = "after",
    inclusive = FALSE
  )
#> # A tibble: 3 × 3
#>   USUBJID AVISITN AVALC
#>   <chr>     <dbl> <chr>
#> 1 1             5 NE   
#> 2 2             2 PD   
#> 3 2             3 PD   

# Select observations from first response to first PD
response %>%
  filter_relative(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    condition = AVALC %in% c("CR", "PR"),
    mode = "first",
    selection = "after",
    inclusive = TRUE,
    keep_no_ref_groups = FALSE
  ) %>%
  filter_relative(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )
#> # A tibble: 7 × 3
#>   USUBJID AVISITN AVALC
#>   <chr>     <dbl> <chr>
#> 1 1             1 PR   
#> 2 1             2 CR   
#> 3 1             3 CR   
#> 4 1             4 SD   
#> 5 1             5 NE   
#> 6 4             2 PR   
#> 7 4             3 PD   
```
