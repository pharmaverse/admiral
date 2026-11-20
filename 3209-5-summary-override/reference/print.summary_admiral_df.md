# Print a Summary of an `admiral_df` Dataset

Print a Summary of an `admiral_df` Dataset

## Usage

``` r
# S3 method for class 'summary_admiral_df'
print(x, ...)
```

## Arguments

- x:

  A `summary_admiral_df` object created by
  [`summary.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)

  Default value

  :   none

- ...:

  Not used

  Default value

  :   none

## Value

`x`, invisibly

## See also

Utilities for printing:
[`print.adam_templates()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.adam_templates.md),
[`print.admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.admiral_df.md),
[`print.duplicates()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.duplicates.md),
[`print.source()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.source.md),
[`print_named_list()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print_named_list.md)

## Examples

``` r
library(tibble)

adsl <- tribble(
  ~STUDYID, ~USUBJID, ~TRT01P,   ~SAFFL,
  "P",      "P-1",    "Placebo", "Y",
  "P",      "P-2",    "Drug A",  "Y",
  "P",      "P-3",    "Drug A",  "N"
)

print(summary(admiral_df(adsl, dataset_keys = exprs(USUBJID), dataset_name = "ADSL")))
#> ── ADSL summary ────────────────────────────────────────────────────────────────
#> Subjects (USUBJID): 3
#> Observations: 3 | Variables: 4
#> ✔ Structure: one record per USUBJID
#> Subjects per treatment arm (TRT01P): Drug A: 2 | Placebo: 1
#> Populations: SAFFL 2
```
