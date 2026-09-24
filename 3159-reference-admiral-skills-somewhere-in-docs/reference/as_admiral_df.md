# Tag a Dataset with the `admiral_df` Class

Adds the `admiral_df` class to a data frame (unless it is already
present), preserving the existing classes such as `tbl_df`/`data.frame`.
This is the tagging primitive that admiral's `derive_*()` and related
dataset-producing functions apply to their output, so that a dataset
built with admiral can be recognized as such downstream (e.g. by a
future dedicated [`summary()`](https://rdrr.io/r/base/summary.html)
method).

## Usage

``` r
as_admiral_df(dataset)
```

## Arguments

- dataset:

  A data frame (or `NULL`)

  Default value

  :   none

## Value

If `dataset` is `NULL`, `NULL` is returned unchanged. Otherwise
`dataset` with `"admiral_df"` prepended to its class attribute.

## See also

Utilities used within Derivation functions:
[`extract_unit()`](https:/pharmaverse.github.io/admiral/3159-reference-admiral-skills-somewhere-in-docs/reference/extract_unit.md),
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/3159-reference-admiral-skills-somewhere-in-docs/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/3159-reference-admiral-skills-somewhere-in-docs/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/3159-reference-admiral-skills-somewhere-in-docs/reference/get_vars_query.md)

## Examples

``` r
library(tibble)

adsl <- tribble(
  ~USUBJID, ~AGE,
  "1",      63,
  "2",      71
)

class(as_admiral_df(adsl))
#> [1] "admiral_df" "tbl_df"     "tbl"        "data.frame"
```
