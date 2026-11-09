# Construct an `admiral_df` Object

Low-level constructor of the `admiral_df` class: adds the class to a
data frame (unless it is already present), preserving the existing
classes such as `tbl_df`/`data.frame`, and stores any metadata supplied.
The arguments are not checked, so this is what admiral's `derive_*()`
and related dataset-producing functions apply to their output. Users
call
[`admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/admiral_df.md),
which validates its arguments before calling this.

## Usage

``` r
new_admiral_df(
  x,
  dataset_keys = NULL,
  dataset_name = NULL,
  dataset_class = NULL
)
```

## Arguments

- x:

  A data frame (or `NULL`)

  Default value

  :   none

- dataset_keys:

  Character vector of key variables, stored in the
  `"admiral_dataset_keys"` attribute

  Default value

  :   `NULL`

- dataset_name:

  Character scalar, stored in the `"admiral_dataset_name"` attribute

  Default value

  :   `NULL`

- dataset_class:

  Character scalar, stored in the `"admiral_dataset_class"` attribute

  Default value

  :   `NULL`

## Value

If `x` is `NULL`, `NULL` is returned unchanged. Otherwise `x` with the
supplied attributes set and `"admiral_df"` prepended to its class
attribute.

## Details

An argument left as `NULL` leaves the corresponding attribute untouched,
so that tagging a dataset which already carries metadata does not remove
it.

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/main/reference/admiral-package.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/extract_duplicate_records.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/format.basket_select.md),
[`infer_class_by_structure()`](https:/pharmaverse.github.io/admiral/main/reference/infer_class_by_structure.md),
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/signal_duplicate_records.md)
