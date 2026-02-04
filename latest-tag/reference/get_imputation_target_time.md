# Get Time Imputation Targets

Get Time Imputation Targets

## Usage

``` r
get_imputation_target_time(time_imputation)
```

## Arguments

- time_imputation:

  The value to impute the time when a timepart is missing.

  A character value is expected, either as a

  - format with hour, min and sec specified as `"hh:mm:ss"`: e.g.
    `"00:00:00"` for the start of the day,

  - or as a keyword: `"first"`,`"last"` to impute to the start/end of a
    day.

  Default value

  :   none

## Value

A list of character vectors. The elements of the list are named "hour",
"minute", "second".

## Details

- For `time_imputation = "first"` `"00"`, `"00"`, `"00"` are returned.

- For `time_imputation = "last"` `"23"`, `"59"`, `"59"` are returned.

- For `time_imputation = "<hh>:<mm>:<ss>"` `"<hh>"`, `"<mm>"`, `"<ss>"`
  are returned.

## See also

[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/impute_dtc_dtm.md)

Utilities used for date imputation:
[`dt_level()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/dt_level.md),
[`dtm_level()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/dtm_level.md),
[`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_imputation_target_date.md),
[`get_partialdatetime()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_partialdatetime.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_imputed_dtc_dt.md),
[`restrict_imputed_dtc_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_imputed_dtc_dtm.md)

## Examples

``` r
# Get imputation target for "first" time
target_first_time <- admiral:::get_imputation_target_time("first")
print(target_first_time)
#> $hour
#> [1] "00"
#> 
#> $minute
#> [1] "00"
#> 
#> $second
#> [1] "00"
#> 

# Get imputation target for "last" time
target_last_time <- admiral:::get_imputation_target_time("last")
print(target_last_time)
#> $hour
#> [1] "23"
#> 
#> $minute
#> [1] "59"
#> 
#> $second
#> [1] "59"
#> 

# Get imputation target for custom time imputation "12-34-56"
target_custom_time <- admiral:::get_imputation_target_time("12-34-56")
print(target_custom_time)
#> $hour
#> [1] "12"
#> 
#> $minute
#> [1] "34"
#> 
#> $second
#> [1] "56"
#> 
```
