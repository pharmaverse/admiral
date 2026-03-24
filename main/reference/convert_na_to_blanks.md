# Convert NAs Into Blank Strings

Turn `NA`s to blank strings .

## Usage

``` r
convert_na_to_blanks(x)

# Default S3 method
convert_na_to_blanks(x)

# S3 method for class 'character'
convert_na_to_blanks(x)

# S3 method for class 'list'
convert_na_to_blanks(x)

# S3 method for class 'data.frame'
convert_na_to_blanks(x)
```

## Arguments

- x:

  Any R object

  Default value

  :   none

## Value

An object of the same class as the input

## Details

The default methods simply returns its input unchanged. The `character`
method turns every instance of `NA_character_` or `NA` into `""` while
preserving *all* attributes. When given a data frame as input the
function keeps all non-character columns as is and applies the just
described logic to `character` all attributes such as labels are
preserved.

## See also

Utilities for Formatting Observations:
[`convert_blanks_to_na()`](https:/pharmaverse.github.io/admiral/main/reference/convert_blanks_to_na.md),
[`yn_to_numeric()`](https:/pharmaverse.github.io/admiral/main/reference/yn_to_numeric.md)

## Examples

``` r
library(tibble)

convert_na_to_blanks(c("a", "b", NA, "d", NA))
#> [1] "a" "b" ""  "d" "" 

df <- tribble(
  ~USUBJID,   ~RFICDTC,
  "1001", "2000-01-01",
  "1002", "2001-01-01",
  "1003",           NA
)
print(df)
#> # A tibble: 3 × 2
#>   USUBJID RFICDTC   
#>   <chr>   <chr>     
#> 1 1001    2000-01-01
#> 2 1002    2001-01-01
#> 3 1003    NA        
convert_na_to_blanks(df)
#> # A tibble: 3 × 2
#>   USUBJID RFICDTC     
#>   <chr>   <chr>       
#> 1 1001    "2000-01-01"
#> 2 1002    "2001-01-01"
#> 3 1003    ""          
```
