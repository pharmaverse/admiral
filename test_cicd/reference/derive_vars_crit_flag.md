# Derive Criterion Flag Variables `CRITy`, `CRITyFL`, and `CRITyFN`

The function derives ADaM compliant criterion flags, e.g., to facilitate
subgroup analyses.

If a criterion flag can't be derived with this function, the derivation
is not ADaM compliant. It helps to ensure that:

- the condition of the criterion depends only on variables of the same
  row,

- the `CRITyFL` is populated with valid values, i.e, either `"Y"` and
  `NA` or `"Y"`, `"N"`, and `NA`,

- the `CRITy` variable is populated correctly, i.e.,

  - set to a constant value within a parameter if `CRITyFL` is populated
    with `"Y"`, `"N"`, and `NA` and

  - set to a constant value within a parameter if the criterion
    condition is fulfilled and to `NA` otherwise if `CRITyFL` is
    populated with `"Y"`, and `NA`

## Usage

``` r
derive_vars_crit_flag(
  dataset,
  crit_nr = 1,
  condition,
  description,
  values_yn = FALSE,
  create_numeric_flag = FALSE
)
```

## Arguments

- dataset:

  Input dataset

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- crit_nr:

  The criterion number, i.e., the `y` in `CRITy`

  Permitted values

  :   a positive integer, e.g. `2` or `5`

  Default value

  :   `1`

- condition:

  Condition for flagging records

  See description of the `values_yn` argument for details on how the
  `CRITyFL` variable is populated.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   none

- description:

  The description of the criterion

  The `CRITy` variable is set to the specified value.

  An expression can be specified to set the value depending on the
  parameter. Please note that the value must be constant within a
  parameter.

  Permitted values

  :   an unquoted expression which evaluates to a character (in
      `dataset`)

  Default value

  :   none

- values_yn:

  Should `"Y"` and `"N"` be used for `CRITyFL`?

  If set to `TRUE`, the `CRITyFL` variable is set to `"Y"` if the
  condition (`condition`) evaluates to `TRUE`, it is set to `"N"` if the
  condition evaluate to `FALSE`, and to `NA` if it evaluates to `NA`.

  Otherwise, the `CRITyFL` variable is set to `"Y"` if the condition
  (`condition`) evaluates to `TRUE`, and to `NA` otherwise.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- create_numeric_flag:

  Create a numeric flag?

  If set to `TRUE`, the `CRITyFN` variable is created. It is set to `1`
  if `CRITyFL == "Y"`, it set to `0` if `CRITyFL == "N"`, and to `NA`
  otherwise.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

## Value

The input dataset with the variables `CRITy`, `CRITyFL`, and optionally
`CRITyFN` added.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md)

## Examples

### Data setup

The following examples use the BDS dataset below as a basis.

    library(tibble, warn.conflicts = FALSE)

    adbds <- tribble(
      ~PARAMCD, ~AVAL,
      "AST",    42,
      "AST",    52,
      "AST",    NA_real_,
      "ALT",    33,
      "ALT",    51
    )

### Creating a simple criterion flag with values `"Y"` and `NA` (`condition`, `description`)

The following call is a simple application of `derive_vars_crit_flag()`
to derive a criterion flag/variable pair in a BDS dataset.

- The new variables are named `CRIT1`/`CRIT1FL` because the argument
  `crit_nr` has not been passed.

- Since the argument `values_yn` has also not been passed and thus is
  set to its default of `FALSE`, `CRIT1FL` is set to `Y` only if
  `condition` evaluates to `TRUE`. For example, in both the first and
  third records, where `condition` is respectively `FALSE` and `NA`, we
  set `CRIT1FL = NA_character_`. The fourth record also exhibits this
  behavior. Also, as per CDISC standards, in this case `CRIT1` is
  populated only for records where `condition` evaluates to `TRUE`.

    derive_vars_crit_flag(
      adbds,
      condition = AVAL > 50,
      description = "Absolute value > 50"
    )
    #> # A tibble: 5 × 4
    #>   PARAMCD  AVAL CRIT1FL CRIT1
    #>   <chr>   <dbl> <chr>   <chr>
    #> 1 AST        42 <NA>    <NA>
    #> 2 AST        52 Y       Absolute value > 50
    #> 3 AST        NA <NA>    <NA>
    #> 4 ALT        33 <NA>    <NA>
    #> 5 ALT        51 Y       Absolute value > 50

The `description` argument also accepts expressions which depend on
other variables in the input dataset. This can be useful to dynamically
populate `CRITx`, for instance in the case below where we improve the
`CRIT1` text because the same flag/variable pair is actually being used
for multiple parameters.

    derive_vars_crit_flag(
      adbds,
      condition = AVAL > 50,
      description = paste(PARAMCD, "> 50"),
    )
    #> # A tibble: 5 × 4
    #>   PARAMCD  AVAL CRIT1FL CRIT1
    #>   <chr>   <dbl> <chr>   <chr>
    #> 1 AST        42 <NA>    <NA>
    #> 2 AST        52 Y       AST > 50
    #> 3 AST        NA <NA>    <NA>
    #> 4 ALT        33 <NA>    <NA>
    #> 5 ALT        51 Y       ALT > 50

### Creating a criterion flag with values `"Y"`, `"N"` and `NA` (`values_yn`)

The next call builds on the previous example by using `value_yn = TRUE`
to distinguish between the cases where `condition` is `FALSE` and those
where it is not evaluable at all.

- As compared to the previous example, for the first record `condition`
  evaluates to `FALSE` and so we set `CRIT1FL = "N"`, whereas for the
  third record, `condition` evaluates to `NA` because `AVAL` is missing
  and so we set `CRIT1FL` to `NA`.

- Note also that because we are using the values `"Y"`, `"N"` and `NA`
  for the flag, as per CDISC standards `CRIT1` is now populated for all
  records rather than just for the `"Y"` records.

    derive_vars_crit_flag(
      adbds,
      condition = AVAL > 50,
      description = paste(PARAMCD, "> 50"),
      values_yn = TRUE
    )
    #> # A tibble: 5 × 4
    #>   PARAMCD  AVAL CRIT1FL CRIT1
    #>   <chr>   <dbl> <chr>   <chr>
    #> 1 AST        42 N       AST > 50
    #> 2 AST        52 Y       AST > 50
    #> 3 AST        NA <NA>    AST > 50
    #> 4 ALT        33 N       ALT > 50
    #> 5 ALT        51 Y       ALT > 50

If the user wishes to set the criterion flag to `"N"` whenever the
condition is not fulfilled, `condition` can be updated using an
`if_else` call, where the third argument determines the behavior when
the condition is not evaluable.

    derive_vars_crit_flag(
      adbds,
      condition = if_else(AVAL > 50, TRUE, FALSE, FALSE),
      description = paste(PARAMCD, "> 50"),
      values_yn = TRUE
    )
    #> # A tibble: 5 × 4
    #>   PARAMCD  AVAL CRIT1FL CRIT1
    #>   <chr>   <dbl> <chr>   <chr>
    #> 1 AST        42 N       AST > 50
    #> 2 AST        52 Y       AST > 50
    #> 3 AST        NA N       AST > 50
    #> 4 ALT        33 N       ALT > 50
    #> 5 ALT        51 Y       ALT > 50

### Specifying the criterion variable/flag number and creating a numeric flag (`crit_nr`, `create_numeric_flag`).

The user can manually specify the criterion variable/flag number to use
to name `CRITy`/`CRITyFL` by passing the `crit_nr` argument. This may be
necessary if, for instance, other criterion flags already exist in the
input dataset.

The user can also choose to create an additional, equivalent numeric
flag `CRITyFN` by setting `create_numeric_flag` to `TRUE`.

    derive_vars_crit_flag(
      adbds,
      condition = AVAL > 50,
      description = paste(PARAMCD, "> 50"),
      values_yn = TRUE,
      crit_nr = 2,
      create_numeric_flag = TRUE
    )
    #> # A tibble: 5 × 5
    #>   PARAMCD  AVAL CRIT2FL CRIT2    CRIT2FN
    #>   <chr>   <dbl> <chr>   <chr>      <int>
    #> 1 AST        42 N       AST > 50       0
    #> 2 AST        52 Y       AST > 50       1
    #> 3 AST        NA <NA>    AST > 50      NA
    #> 4 ALT        33 N       ALT > 50       0
    #> 5 ALT        51 Y       ALT > 50       1
