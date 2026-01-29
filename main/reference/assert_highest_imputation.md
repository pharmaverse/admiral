# Assert Highest Imputation Validity

This function checks the validity and requirements for the
`highest_imputation` argument. It ensures that necessary conditions for
`date_imputation`, `min_dates`, and `max_dates` are met when
`highest_imputation` is set to `"Y"`.

## Usage

``` r
assert_highest_imputation(
  highest_imputation,
  highest_imputation_values,
  date_imputation = NULL,
  max_dates,
  min_dates
)
```

## Arguments

- highest_imputation:

  A character scalar indicating the highest level of imputation.

  Default value

  :   none

- highest_imputation_values:

  A character vector of valid values for `highest_imputation`.

  Default value

  :   none

- date_imputation:

  Optional character scalar specifying the imputation method for dates.

  Default value

  :   `NULL`

- max_dates:

  Optional vector specifying maximum dates for imputation.

  Default value

  :   none

- min_dates:

  Optional vector specifying minimum dates for imputation.

  Default value

  :   none

## Value

Returns `NULL` invisibly if assertions pass.

## Details

- If `highest_imputation` is "Y", either `min_dates` or `max_dates` must
  be specified.

- If `highest_imputation` is "Y" and `date_imputation` is "first",
  `min_dates` must be specified.

- If `highest_imputation` is "Y" and `date_imputation` is "last",
  `max_dates` must be specified.
