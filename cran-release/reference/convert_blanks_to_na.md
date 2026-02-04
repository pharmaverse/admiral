# Convert Blank Strings Into NAs

Turn SAS blank strings into proper R `NA`s.

## Usage

``` r
convert_blanks_to_na(x)

# Default S3 method
convert_blanks_to_na(x)

# S3 method for class 'character'
convert_blanks_to_na(x)

# S3 method for class 'list'
convert_blanks_to_na(x)

# S3 method for class 'data.frame'
convert_blanks_to_na(x)
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
method turns every instance of `""` into `NA_character_` while
preserving *all* attributes. When given a data frame as input the
function keeps all non-character columns as is and applies the just
described logic to `character` columns. Once again all attributes such
as labels are preserved.

## See also

Utilities for Formatting Observations:
[`convert_na_to_blanks()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/convert_na_to_blanks.md),
[`yn_to_numeric()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/yn_to_numeric.md)

## Examples

``` r
library(tibble)

convert_blanks_to_na(c("a", "b", "", "d", ""))
#> [1] "a" "b" NA  "d" NA 

df <- tribble(
  ~USUBJID,   ~RFICDTC,
  "1001", "2000-01-01",
  "1002", "2001-01-01",
  "1003",           ""
)
print(df)
#> # A tibble: 3 × 2
#>   USUBJID RFICDTC     
#>   <chr>   <chr>       
#> 1 1001    "2000-01-01"
#> 2 1002    "2001-01-01"
#> 3 1003    ""          
convert_blanks_to_na(df)
#> # A tibble: 3 × 2
#>   USUBJID RFICDTC   
#>   <chr>   <chr>     
#> 1 1001    2000-01-01
#> 2 1002    2001-01-01
#> 3 1003    NA        
```
