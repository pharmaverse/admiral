# Format Imputed Date/Datetime

Formats imputed date/datetime components into a string representation.

## Usage

``` r
format_imputed_dtc(imputed)
```

## Arguments

- imputed:

  A list of imputed date/time components.

  Default value

  :   none

## Value

A character vector of formatted date/datetime strings.

## Details

The function formats the imputed components into `"YYYY-MM-DD"` for
dates and `"YYYY-MM-DDThh:mm:ss"` for datetimes. It replaces any string
containing `"x"` with `NA`.

## Examples

``` r
# Format imputed datetime components
imputed_datetime <- list(
  year = "2020", month = "01", day = "01",
  hour = "12", minute = "00", second = "00"
)
formatted_datetime <- admiral:::format_imputed_dtc(imputed_datetime)
print(formatted_datetime)
#> [1] "2020-01-01T12:00:00"

# Format imputed date components
imputed_date <- list(year = "2020", month = "01", day = "01")
formatted_date <- admiral:::format_imputed_dtc(imputed_date)
print(formatted_date)
#> [1] "2020-01-01"

# Handle imputed datetime with missing parts (contains 'x')
# Expected: NA because 'x' is an undefined component
imputed_partial_datetime <- list(
  year = "2020", month = "xx", day = "01",
  hour = "12", minute = "00", second = "00"
)
formatted_partial_datetime <- admiral:::format_imputed_dtc(imputed_partial_datetime)
print(formatted_partial_datetime)
#> [1] NA

# Handle imputed date with missing parts (contains 'x')
# Expected: NA because 'x' is an undefined component
imputed_partial_date <- list(year = "2020", month = "xx", day = "01")
formatted_partial_date <- admiral:::format_imputed_dtc(imputed_partial_date)
print(formatted_partial_date)
#> [1] NA
```
