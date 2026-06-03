# Get Date Imputation Targets

Get Date Imputation Targets

## Usage

``` r
get_imputation_target_date(date_imputation, month)
```

## Arguments

- date_imputation:

  The value to impute the day/month when a datepart is missing.

  A character value is expected, as a:

  - format with month and day specified as `"mm-dd"`: e.g. `"06-15"` for
    the 15th of June,

  - format with day specified as `"dd"`: e.g. `"15"` for the 15th day of
    a month

  - or as a keyword: `"first"`, `"mid"`, `"last"` to impute to the
    first/mid/last day/month.

  Default value

  :   none

- month:

  Month component of the partial date

  Default value

  :   none

## Value

A list of character vectors. The elements of the list are named "year",
"month", "day".

## Details

- For `date_imputation = "first"`, `"0000"`, `"01"`, `"01"` are
  returned.

- For `date_imputation = "mid"`, `"xxxx"`, `"06"`, `"30"` if `month` is
  `NA`. otherwise `"15"` returned.

- For `date_imputation = "last"`, `"9999"`, `"12"`, `"28"` are returned.

- For `date_imputation = "<mm>-<dd>"`, `"xxxx"`, `"<mm>"`, `"<dd>"` are
  returned.

- For `date_imputation = "<dd>"`, `"xxxx"`, `"xx"`, `"<dd>"` are
  returned.

`"xxxx"` indicates that the component is undefined. If an undefined
component occurs in the imputed `--DTC` value, the imputed `--DTC` value
is set to `NA_character_` in the imputation functions.

## See also

[`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/impute_dtc_dtm.md),
[`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/impute_dtc_dt.md)

Utilities used for date imputation:
[`dt_level()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/dt_level.md),
[`dtm_level()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/dtm_level.md),
[`get_imputation_target_time()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/get_imputation_target_time.md),
[`get_partialdatetime()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/get_partialdatetime.md),
[`restrict_imputed_dtc_dt()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/restrict_imputed_dtc_dt.md),
[`restrict_imputed_dtc_dtm()`](https:/pharmaverse.github.io/admiral/3014_derive_vars_dt_user_defined_date_imp_with_d/reference/restrict_imputed_dtc_dtm.md)

## Examples

``` r
# Get imputation target for "first"
admiral:::get_imputation_target_date("first", month = NA)
#> $year
#> [1] "0000"
#> 
#> $month
#> [1] "01"
#> 
#> $day
#> [1] "01"
#> 

# Get imputation target for "mid" with specified month
admiral:::get_imputation_target_date("mid", month = "03")
#> $year
#> [1] "xxxx"
#> 
#> $month
#> [1] "06"
#> 
#> $day
#> [1] "15"
#> 

# Get imputation target for "mid" with NA month
admiral:::get_imputation_target_date("mid", month = NA)
#> $year
#> [1] "xxxx"
#> 
#> $month
#> [1] "06"
#> 
#> $day
#> [1] "30"
#> 

# Get imputation target for "last"
admiral:::get_imputation_target_date("last", month = NA)
#> $year
#> [1] "9999"
#> 
#> $month
#> [1] "12"
#> 
#> $day
#> [1] "28"
#> 

# Get imputation target for custom date imputation "06-15"
admiral:::get_imputation_target_date("06-15", month = NA)
#> $year
#> [1] "xxxx"
#> 
#> $month
#> [1] "06"
#> 
#> $day
#> [1] "15"
#> 

# Get imputation target for custom date imputation "11"
admiral:::get_imputation_target_date("11", month = NA)
#> $year
#> [1] "xxxx"
#> 
#> $month
#> [1] "xx"
#> 
#> $day
#> [1] "11"
#> 
```
