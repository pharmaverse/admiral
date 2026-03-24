# Get One to Many Values that Led to a Prior Error

Get One to Many Values that Led to a Prior Error

## Usage

``` r
get_one_to_many_dataset()
```

## Value

A `data.frame` or `NULL`

## Details

If
[`assert_one_to_one()`](https://pharmaverse.github.io/admiraldev/reference/assert_one_to_one.html)
detects an issue, the one to many values are stored in a dataset. This
dataset can be retrieved by `get_one_to_many_dataset()`.

Note that the function always returns the one to many values from the
last error that has been thrown in the current R session. Thus, after
restarting the R sessions `get_one_to_many_dataset()` will return `NULL`
and after a second error has been thrown, the dataset of the first error
can no longer be accessed (unless it has been saved in a variable).

## See also

Utilities for Dataset Checking:
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/get_duplicates_dataset.md),
[`get_many_to_one_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/get_many_to_one_dataset.md)

## Examples

``` r
library(admiraldev, warn.conflicts = FALSE)
data(admiral_adsl)

try(
  assert_one_to_one(admiral_adsl, exprs(STUDYID), exprs(SITEID))
)
#> Error in eval(expr, envir) : 
#>   For some values of "STUDYID" there is more than one value of "SITEID"
#> ℹ Call `admiral::get_one_to_many_dataset()` to get all one-to-many values.

get_one_to_many_dataset()
#> # A tibble: 17 × 2
#> # Groups:   STUDYID [1]
#>    STUDYID      SITEID
#>    <chr>        <chr> 
#>  1 CDISCPILOT01 701   
#>  2 CDISCPILOT01 702   
#>  3 CDISCPILOT01 703   
#>  4 CDISCPILOT01 704   
#>  5 CDISCPILOT01 705   
#>  6 CDISCPILOT01 706   
#>  7 CDISCPILOT01 707   
#>  8 CDISCPILOT01 708   
#>  9 CDISCPILOT01 709   
#> 10 CDISCPILOT01 710   
#> 11 CDISCPILOT01 711   
#> 12 CDISCPILOT01 713   
#> 13 CDISCPILOT01 714   
#> 14 CDISCPILOT01 715   
#> 15 CDISCPILOT01 716   
#> 16 CDISCPILOT01 717   
#> 17 CDISCPILOT01 718   
```
