# Check if a Partial Date/Time is a Datetime

This function determines whether a given partial date/time structure
represents a datetime or just a date.

## Usage

``` r
is_partial_datetime(partial)
```

## Arguments

- partial:

  A named list containing date or datetime components.

  Default value

  :   none

## Value

A logical value. TRUE if the partial represents a datetime, FALSE if it
represents a date only.

## Details

The function checks for the presence of all date components (year,
month, day) and all time components (hour, minute, second) in the input
list. If all components are present, it's considered a datetime. If only
date components are present, it's considered a date. Any other
combination will result in an error.

## Examples

``` r
# Datetime example
partial_datetime <- list(
  year = "2023", month = "05", day = "15",
  hour = "14", minute = "30", second = "00"
)
admiral:::is_partial_datetime(partial_datetime) # Returns TRUE
#> [1] TRUE

# Date example
partial_date <- list(year = "2023", month = "05", day = "15")
admiral:::is_partial_datetime(partial_date) # Returns FALSE
#> [1] FALSE

# Invalid example
if (FALSE) { # \dontrun{
partial_invalid <- list(year = "2023", month = "05", hour = "14")
admiral:::is_partial_datetime(partial_invalid) # Throws an error
} # }
```
