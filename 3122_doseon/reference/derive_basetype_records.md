# Derive Basetype Variable

Baseline Type `BASETYPE` is needed when there is more than one
definition of baseline for a given Analysis Parameter `PARAM` in the
same dataset. For a given parameter, if Baseline Value `BASE` or `BASEC`
are derived and there is more than one definition of baseline, then
`BASETYPE` must be non-null on all records of any type for that
parameter where either `BASE` or `BASEC` are also non-null. Each value
of `BASETYPE` refers to a definition of baseline that characterizes the
value of `BASE` on that row. Please see section 4.2.1.6 of the ADaM
Implementation Guide, version 1.3 for further background.

## Usage

``` r
derive_basetype_records(dataset, basetypes)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `basetypes` argument are expected to be
  in the dataset.

  Default value

  :   none

- basetypes:

  A *named* list of expressions created using the
  [`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  function

  The names corresponds to the values of the newly created `BASETYPE`
  variables and the expressions are used to subset the input dataset.

  Default value

  :   none

## Value

The input dataset with variable `BASETYPE` added

## Details

Adds the `BASETYPE` variable to a dataset and duplicates records based
upon the provided conditions.

For each element of `basetypes` the input dataset is subset based upon
the provided expression and the `BASETYPE` variable is set to the name
of the expression. Then, all subsets are stacked. Records which do not
match any condition are kept and `BASETYPE` is set to `NA`.

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/3122_doseon/reference/derive_summary_records.md)

## Examples

### Add records for different baseline types (`basetypes`)

The `basetypes` argument is a named list of expressions where each name
becomes a value of `BASETYPE` and each expression defines which records
receive that value. A record can match multiple expressions and will be
duplicated once for each matching `BASETYPE`. In this example, records
for subject `P01` show the duplication across both baseline types.

Records that do not match any condition in `basetypes` are kept in the
output dataset with `BASETYPE` set to `NA`. In this example, `SCREENING`
records do not match any of the `basetypes` conditions and are therefore
retained with `BASETYPE = NA`.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)

    bds <- tribble(
      ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
      "P01",    "SCREENING",    "PARAM01",     1,  10.2,
      "P01",    "RUN-IN",       "PARAM01",     2,  10.0,
      "P01",    "RUN-IN",       "PARAM01",     3,   9.8,
      "P01",    "DOUBLE-BLIND", "PARAM01",     4,   9.2,
      "P01",    "DOUBLE-BLIND", "PARAM01",     5,  10.1,
      "P02",    "SCREENING",    "PARAM01",     1,  12.2,
      "P02",    "RUN-IN",       "PARAM01",     2,  12.1,
      "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.2
    )

    derive_basetype_records(
      dataset = bds,
      basetypes = exprs(
        "RUN-IN" = EPOCH %in% c("RUN-IN", "DOUBLE-BLIND"),
        "DOUBLE-BLIND" = EPOCH == "DOUBLE-BLIND"
      )
    )
    #> # A tibble: 11 × 6
    #>    USUBJID EPOCH        PARAMCD  ASEQ  AVAL BASETYPE
    #>    <chr>   <chr>        <chr>   <dbl> <dbl> <chr>
    #>  1 P01     SCREENING    PARAM01     1  10.2 <NA>
    #>  2 P02     SCREENING    PARAM01     1  12.2 <NA>
    #>  3 P01     RUN-IN       PARAM01     2  10   RUN-IN
    #>  4 P01     RUN-IN       PARAM01     3   9.8 RUN-IN
    #>  5 P01     DOUBLE-BLIND PARAM01     4   9.2 RUN-IN
    #>  6 P01     DOUBLE-BLIND PARAM01     5  10.1 RUN-IN
    #>  7 P02     RUN-IN       PARAM01     2  12.1 RUN-IN
    #>  8 P02     DOUBLE-BLIND PARAM01     3  10.2 RUN-IN
    #>  9 P01     DOUBLE-BLIND PARAM01     4   9.2 DOUBLE-BLIND
    #> 10 P01     DOUBLE-BLIND PARAM01     5  10.1 DOUBLE-BLIND
    #> 11 P02     DOUBLE-BLIND PARAM01     3  10.2 DOUBLE-BLIND

### Include all records for multiple baseline type derivations (`basetypes = TRUE`)

When all parameter records need to be included for multiple baseline
type derivations (such as `"LAST"` and `"WORST"`), set each expression
in `basetypes` to `TRUE`. This duplicates every record once for each
named baseline type.

    bds <- tribble(
      ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
      "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
      "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
      "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
      "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1
    )

    derive_basetype_records(
      dataset = bds,
      basetypes = exprs(
        "LAST" = TRUE,
        "WORST" = TRUE
      )
    )
    #> # A tibble: 8 × 6
    #>   USUBJID EPOCH        PARAMCD  ASEQ  AVAL BASETYPE
    #>   <chr>   <chr>        <chr>   <dbl> <dbl> <chr>
    #> 1 P01     RUN-IN       PARAM01     1  10   LAST
    #> 2 P01     RUN-IN       PARAM01     2   9.8 LAST
    #> 3 P01     DOUBLE-BLIND PARAM01     3   9.2 LAST
    #> 4 P01     DOUBLE-BLIND PARAM01     4  10.1 LAST
    #> 5 P01     RUN-IN       PARAM01     1  10   WORST
    #> 6 P01     RUN-IN       PARAM01     2   9.8 WORST
    #> 7 P01     DOUBLE-BLIND PARAM01     3   9.2 WORST
    #> 8 P01     DOUBLE-BLIND PARAM01     4  10.1 WORST   
