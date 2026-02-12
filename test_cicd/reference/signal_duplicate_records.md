# Signal Duplicate Records

Signal Duplicate Records

## Usage

``` r
signal_duplicate_records(
  dataset,
  by_vars,
  msg = paste("Dataset contains duplicate records", "with respect to",
    "{.var {replace_values_by_names(by_vars)}}"),
  cnd_type = "error",
  class = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Defines groups of records in which to look for duplicates.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- msg:

  The condition message

  Default value

  :   `paste("Dataset contains duplicate records", "with respect to", "{.var {replace_values_by_names(by_vars)}}")`

- cnd_type:

  Type of condition to signal when detecting duplicate records.

  Permitted values

  :   `"message"`, `"warning"`, or `"error"`

  Default value

  :   `"error"`

- class:

  Class of the condition

  The specified classes are added to the classes of the condition.
  `c("duplicate_records", "assert-admiral")` is always added.

  Default value

  :   `NULL`

## Value

No return value, called for side effects

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/test_cicd/reference/admiral-package.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/extract_duplicate_records.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/format.basket_select.md)

## Examples

``` r
data(admiral_adsl)

# Duplicate the first record
adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)

signal_duplicate_records(adsl, exprs(USUBJID), cnd_type = "message")
#> Dataset contains duplicate records with respect to `USUBJID`
#> ℹ Run `admiral::get_duplicates_dataset()` to access the duplicate records
```
