# Impute Numeric Values from Character Results with Qualifiers

Converts character results to numeric values, with an optional shift for
strict inequality qualifiers.

## Usage

``` r
compute_qual_imputation(
  character_value,
  imputation_type = "strip",
  factor = NULL,
  offset = 0
)
```

## Arguments

- character_value:

  Character vector of numeric results, optionally prefixed with a
  qualifier, for example, `c("<40", ">=40.1", "3")`.

  Default value

  :   none

- imputation_type:

  How to convert the character results:

  - `"strip"`: Remove qualifiers and use `offset` as the additive offset
    for strict inequalities.

  - `"shift"`: Remove qualifiers and use the decimal precision of each
    original character result as the additive offset. See Details.

  - **\[deprecated\]** `1` - equivalent to `"strip"`.

  - **\[deprecated\]** `2` - equivalent to `"shift"`.

  Permitted values

  :   `"strip"`, `"shift"`

  Default value

  :   `"strip"`

- factor:

  **\[deprecated\]** Deprecated alias for `offset`.

  Default value

  :   `NULL`

- offset:

  Finite, nonnegative numeric scalar used when
  `imputation_type = "strip"`. It is added for results with `>` and
  subtracted for results with `<`. It is ignored when
  `imputation_type = "shift"`.

  Default value

  :   `0`

## Value

A numeric vector of imputed results.

## Details

Qualified results must match `^[<>]=?[0-9]+[.]?[0-9]*$`, for example,
`"<40"`, `">=40.1"`, or `">40."`. Unqualified decimal numbers,
optionally signed, are also accepted, including `"-3"`, `"+3"`, and
`".5"`. Whitespace, scientific notation, signed qualified values, and
malformed qualifiers are not accepted. Invalid strings and missing
values return `NA_real_` without a warning. Results keep the input order
and length.

For a strict `>` qualifier, add the offset; for a strict `<` qualifier,
subtract it. Results with `>=` or `<=`, or without an inequality
qualifier, are converted without a shift.

With `imputation_type = "shift"`, the offset is \\10^{-d}\\, where \\d\\
is the number of decimal places in the original character result,
including trailing zeros. For integer results, the offset is 1. For
example, `"<40"` becomes 39, `">40.1"` becomes 40.2, and `">40.10"`
becomes 40.11.

Numeric methods `1` and `2` are deprecated aliases for `"strip"` and
`"shift"`, respectively. Both still work and produce a deprecation
message.

Choose an imputation rule that matches the analysis specification. A
shift can cross zero, so results can be negative.

## See also

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/main/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qtc.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/main/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/main/reference/transform_range.md)

## Examples

``` r
values <- c("<40", ">40.1", "<=40", ">=40.1", "3")
compute_qual_imputation(values)
#> [1] 40.0 40.1 40.0 40.1  3.0
compute_qual_imputation(values, offset = 0.5)
#> [1] 39.5 40.6 40.0 40.1  3.0
compute_qual_imputation(values, imputation_type = "shift")
#> [1] 39.0 40.2 40.0 40.1  3.0

# Trailing zeros retain the reported precision
compute_qual_imputation(c(">40.1", ">40.10"), imputation_type = "shift")
#> [1] 40.20 40.11

# Malformed qualifiers produce missing values without warnings
compute_qual_imputation(c("0>1", "-<0.9>", ">100.00", "<1"), "shift")
#> [1]     NA     NA 100.01   0.00
```
