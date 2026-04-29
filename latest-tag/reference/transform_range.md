# Transform Range

Transforms results from the source range to the target range. For
example, for transforming source values 1, 2, 3, 4, 5 to 0, 25, 50, 75,
100.

## Usage

``` r
transform_range(
  source,
  source_range,
  target_range,
  flip_direction = FALSE,
  outside_range = "NA"
)
```

## Arguments

- source:

  A vector of values to be transformed

  A numeric vector is expected.

  Default value

  :   none

- source_range:

  The permitted source range

  A numeric vector containing two elements is expected, representing the
  lower and upper bounds of the permitted source range.

  Default value

  :   none

- target_range:

  The target range

  A numeric vector containing two elements is expected, representing the
  lower and upper bounds of the target range.

  Default value

  :   none

- flip_direction:

  Flip direction of the range?

  The transformed values will be reversed within the target range, e.g.
  within the range 0 to 100, 25 would be reversed to 75.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- outside_range:

  Handling of values outside the source range

  Values outside the source range (`source_range`) are transformed to
  `NA`.

  If `"warning"` or `"error"` is specified, a warning or error is issued
  if `source` includes any values outside the source range.

  Permitted values

  :   `"NA"`, `"warning"`, `"error"`

  Default value

  :   `"NA"`

## Value

The source linearly transformed to the target range

## Details

Returns the values of `source` linearly transformed from the source
range (`source_range`) to the target range (`target_range`). Values
outside the source range are set to `NA`.

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
[`compute_scale()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_scale.md)

## Examples

``` r
transform_range(
  source = c(1, 4, 3, 6, 5),
  source_range = c(1, 5),
  target_range = c(0, 100)
)
#> [1]   0  75  50  NA 100

transform_range(
  source = c(1, 4, 3, 6, 5),
  source_range = c(1, 5),
  target_range = c(0, 100),
  flip_direction = TRUE
)
#> [1] 100  25  50  NA   0
```
