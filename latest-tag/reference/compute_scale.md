# Compute Scale Parameters

Computes the average of a set of source values and transforms the result
from the source range to the target range. For example, for calculating
the average of a set of questionnaire response scores and re-coding the
average response to obtain a subscale score.

## Usage

``` r
compute_scale(
  source,
  source_range = NULL,
  target_range = NULL,
  flip_direction = FALSE,
  min_n = 1
)
```

## Arguments

- source:

  A vector of values to be scaled

  A numeric vector is expected.

  Default value

  :   none

- source_range:

  The permitted source range

  A numeric vector containing two elements is expected, representing the
  lower and upper bounds of the permitted source range. Alternatively,
  if no argument is specified for `source_range` and `target_range`, no
  transformation will be performed.

  Default value

  :   `NULL`

- target_range:

  The target range

  A numeric vector containing two elements is expected, representing the
  lower and upper bounds of the target range. Alternatively, if no
  argument is specified for `source_range` and `target_range`, no
  transformation will be performed.

  Default value

  :   `NULL`

- flip_direction:

  Flip direction of the scale?

  The transformed values will be reversed within the target range, e.g.
  within the range 0 to 100, 25 would be reversed to 75.

  This argument will be ignored if `source_range` and `target_range`
  aren't specified.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- min_n:

  Minimum number of values for computation

  The minimum number of non-missing values in source for the computation
  to be carried out. If the number of non-missing values is below
  `min_n`, the result will be set to missing, i.e. `NA`.

  A positive integer is expected.

  Default value

  :   `1`

## Value

The average of source transformed to the target range or `NA` if source
doesn't contain `min_n` values.

## Details

Returns a numeric value. If source contains less than `min_n` values,
the result is set to `NA`. If `source_range` and `target_range` aren't
specified, the mean will be computed without any transformation being
performed.

## See also

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_rr.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/transform_range.md)

## Examples

``` r
compute_scale(
  source = c(1, 4, 3, 5),
  source_range = c(1, 5),
  target_range = c(0, 100),
  flip_direction = TRUE,
  min_n = 3
)
#> [1] 43.75
```
