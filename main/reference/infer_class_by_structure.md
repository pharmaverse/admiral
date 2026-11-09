# Infer the ADaM Dataset Class of Subject or Device Data

Applies the precedence described in
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md)
for one class of data, subject or device.

## Usage

``` r
infer_class_by_structure(dataset, id_vars, classes)
```

## Arguments

- dataset:

  A data frame

  Default value

  :   none

- id_vars:

  Character vector of the variables identifying the unit of the class
  (the subject keys, or `STUDYID` and `SPDEVID`)

  Default value

  :   none

- classes:

  Named character vector of the class to return for each structure:
  `tte`, `bds`, `occds`, and `level` (one record per unit)

  Default value

  :   none

## Value

A character scalar: one of `classes`, `"other"`, or `NA`.

## See also

Other internal:
[`admiral-package`](https:/pharmaverse.github.io/admiral/main/reference/admiral-package.md),
[`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/extract_duplicate_records.md),
[`format.basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/format.basket_select.md),
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/main/reference/infer_dataset_class.md),
[`new_admiral_df()`](https:/pharmaverse.github.io/admiral/main/reference/new_admiral_df.md),
[`signal_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/signal_duplicate_records.md)
