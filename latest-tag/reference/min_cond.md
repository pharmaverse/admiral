# Minimum Value on a Subset

The function derives the minimum value of a vector/column on a subset of
entries/observations.

## Usage

``` r
min_cond(var, cond)
```

## Arguments

- var:

  A vector

  Default value

  :   none

- cond:

  A condition

  Default value

  :   none

## See also

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/count_vals.md),
[`filter_exist()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_exist.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_extreme.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_joined.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_not_exist.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/max_cond.md)

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
)

# In oncology setting, when needing to check the first time a patient had
# a Complete Response (CR) to compare to see if any Partial Response (PR)
# occurred after this add variable indicating if PR occurred after CR
group_by(data, USUBJID) %>% mutate(
  first_cr_vis = min_cond(var = AVISITN, cond = AVALC == "CR"),
  last_pr_vis = max_cond(var = AVISITN, cond = AVALC == "PR"),
  pr_after_cr = last_pr_vis > first_cr_vis
)
#> # A tibble: 8 × 6
#> # Groups:   USUBJID [2]
#>   USUBJID AVISITN AVALC first_cr_vis last_pr_vis pr_after_cr
#>   <chr>     <dbl> <chr>        <dbl>       <dbl> <lgl>      
#> 1 1             1 PR               2           1 FALSE      
#> 2 1             2 CR               2           1 FALSE      
#> 3 1             3 NE               2           1 FALSE      
#> 4 1             4 CR               2           1 FALSE      
#> 5 1             5 NE               2           1 FALSE      
#> 6 2             1 CR               1           2 TRUE       
#> 7 2             2 PR               1           2 TRUE       
#> 8 2             3 CR               1           2 TRUE       
```
