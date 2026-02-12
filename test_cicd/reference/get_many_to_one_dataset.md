# Get Many to One Values that Led to a Prior Error

Get Many to One Values that Led to a Prior Error

## Usage

``` r
get_many_to_one_dataset()
```

## Value

A `data.frame` or `NULL`

## Details

If
[`assert_one_to_one()`](https://pharmaverse.github.io/admiraldev/reference/assert_one_to_one.html)
detects an issue, the many to one values are stored in a dataset. This
dataset can be retrieved by `get_many_to_one_dataset()`.

Note that the function always returns the many to one values from the
last error that has been thrown in the current R session. Thus, after
restarting the R sessions `get_many_to_one_dataset()` will return `NULL`
and after a second error has been thrown, the dataset of the first error
can no longer be accessed (unless it has been saved in a variable).

## See also

Utilities for Dataset Checking:
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_duplicates_dataset.md),
[`get_one_to_many_dataset()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_one_to_many_dataset.md)

## Examples

``` r
library(admiraldev, warn.conflicts = FALSE)
data(admiral_adsl)

try(
  assert_one_to_one(admiral_adsl, exprs(SITEID), exprs(STUDYID))
)
#> Error in eval(expr, envir) : 
#>   There is more than one value of "SITEID" for some values of "STUDYID"
#> ℹ Call `admiral::get_many_to_one_dataset()` to get all many-to-one values.

get_many_to_one_dataset()
#> # A tibble: 17 × 2
#> # Groups:   STUDYID [1]
#>    SITEID STUDYID     
#>    <chr>  <chr>       
#>  1 701    CDISCPILOT01
#>  2 702    CDISCPILOT01
#>  3 703    CDISCPILOT01
#>  4 704    CDISCPILOT01
#>  5 705    CDISCPILOT01
#>  6 706    CDISCPILOT01
#>  7 707    CDISCPILOT01
#>  8 708    CDISCPILOT01
#>  9 709    CDISCPILOT01
#> 10 710    CDISCPILOT01
#> 11 711    CDISCPILOT01
#> 12 713    CDISCPILOT01
#> 13 714    CDISCPILOT01
#> 14 715    CDISCPILOT01
#> 15 716    CDISCPILOT01
#> 16 717    CDISCPILOT01
#> 17 718    CDISCPILOT01
```
