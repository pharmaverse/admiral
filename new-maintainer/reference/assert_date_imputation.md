# Assert `date_imputation`

Applies assertions on the `date_imputation` argument to reduce
cyclomatic complexity

## Usage

``` r
assert_date_imputation(date_imputation, highest_imputation)
```

## Arguments

- date_imputation:

  The value to impute the day/month when a datepart is missing.

  Default value

  :   none

- highest_imputation:

  Highest imputation level

  Default value

  :   none

## Value

asserted `date_imputation`

## Details

Asserts that date_imputation is a scalar. Asserts that the values in
`date_imputation` are permitted. The permitted values in
`date_imputation` vary by `highest_imputation`
