# Count Number of Observations Where a Variable Equals a Value

Count number of observations where a variable equals a value.

## Usage

``` r
count_vals(var, val)
```

## Arguments

- var:

  A vector

  Default value

  :   none

- val:

  A value

  Default value

  :   none

## See also

Utilities for Filtering Observations:
[`filter_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_exist.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/main/reference/filter_extreme.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_not_exist.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/main/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(admiral)
data <- tribble(
  ~USUBJID, ~AVISITN, ~AVALC,
  "1",      1,        "PR",
  "1",      2,        "CR",
  "1",      3,        "NE",
  "1",      4,        "CR",
  "1",      5,        "NE",
  "2",      1,        "CR",
  "2",      2,        "PR",
  "2",      3,        "CR",
  "3",      1,        "CR",
  "4",      1,        "CR",
  "4",      2,        "NE",
  "4",      3,        "NE",
  "4",      4,        "CR",
  "4",      5,        "PR"
)

# add variable providing the number of NEs for each subject
group_by(data, USUBJID) %>%
  mutate(nr_nes = count_vals(var = AVALC, val = "NE"))
#> # A tibble: 14 × 4
#> # Groups:   USUBJID [4]
#>    USUBJID AVISITN AVALC nr_nes
#>    <chr>     <dbl> <chr>  <int>
#>  1 1             1 PR         2
#>  2 1             2 CR         2
#>  3 1             3 NE         2
#>  4 1             4 CR         2
#>  5 1             5 NE         2
#>  6 2             1 CR         0
#>  7 2             2 PR         0
#>  8 2             3 CR         0
#>  9 3             1 CR         0
#> 10 4             1 CR         2
#> 11 4             2 NE         2
#> 12 4             3 NE         2
#> 13 4             4 CR         2
#> 14 4             5 PR         2
```
