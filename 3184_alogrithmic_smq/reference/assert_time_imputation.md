# Assert `time_imputation`

Applies assertions on the `time_imputation` argument

## Usage

``` r
assert_time_imputation(time_imputation, highest_imputation)
```

## Arguments

- time_imputation:

  The value to impute time when missing

  Default value

  :   none

- highest_imputation:

  Highest imputation level

  Default value

  :   none

## Value

asserted `time_imputation`

## Examples

``` r
# Assert valid 'first' time imputation
time_imp_first <- admiral:::assert_time_imputation("first", "Y")
print(time_imp_first)
#> [1] "first"

# Assert valid 'last' time imputation
time_imp_last <- admiral:::assert_time_imputation("last", "Y")
print(time_imp_last)
#> [1] "last"

# Assert valid custom time imputation "12:34:56"
time_imp_custom <- admiral:::assert_time_imputation("12:34:56", "Y")
print(time_imp_custom)
#> [1] "12:34:56"
```
