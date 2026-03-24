# Derive Categorization Variables Like `AVALCATy` and `AVALCAyN`

Derive Categorization Variables Like `AVALCATy` and `AVALCAyN`

## Usage

``` r
derive_vars_cat(dataset, definition, by_vars = NULL)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `definition` arguments
  are expected to be in the dataset.

  Default value

  :   none

- definition:

  List of expressions created by
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md).
  Must be in rectangular format and specified using the same syntax as
  when creating a `tibble` using the
  [`tribble()`](https://tibble.tidyverse.org/reference/tribble.html)
  function. The `definition` object will be converted to a `tibble`
  using
  [`tribble()`](https://tibble.tidyverse.org/reference/tribble.html)
  inside this function.

  Must contain:

  - the column `condition` which will be converted to a logical
    expression and will be used on the `dataset` input.

  - at least one additional column with the new column name and the
    category value(s) used by the logical expression.

  - the column specified in `by_vars` (if `by_vars` is specified)

  e.g. if `by_vars` is not specified:

      exprs(~condition,   ~AVALCAT1, ~AVALCA1N,
            AVAL >= 140, ">=140 cm",         1,
            AVAL < 140,   "<140 cm",         2)

  e.g. if `by_vars` is specified as `exprs(VSTEST)`:

      exprs(~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
            "Height", AVAL >= 140, ">=140 cm",         1,
            "Height",  AVAL < 140,  "<140 cm",         2)

  Default value

  :   none

- by_vars:

  list of expressions with one element. `NULL` by default. Allows for
  specifying by groups, e.g. `exprs(PARAMCD)`. Variable must be present
  in both `dataset` and `definition`. The conditions in `definition` are
  applied only to those records that match `by_vars`. The categorization
  variables are set to `NA` for records not matching any of the by
  groups in `definition`.

  Default value

  :   `NULL`

## Value

The input dataset with the new variables defined in `definition` added

## Details

If conditions are overlapping, the row order of `definitions` must be
carefully considered. The **first** match will determine the category.
i.e. if

`AVAL = 155`

and the `definition` is:

    definition <- exprs(
      ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
      "Height",  AVAL > 170,  ">170 cm",         1,
      "Height", AVAL <= 170, "<=170 cm",         2,
      "Height", AVAL <= 160, "<=160 cm",         3
    )

then `AVALCAT1` will be `"<=170 cm"`, as this is the first match for
`AVAL`. If you specify:

    definition <- exprs(
      ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
      "Height", AVAL <= 160, "<=160 cm",         3,
      "Height", AVAL <= 170, "<=170 cm",         2,
      "Height",  AVAL > 170,  ">170 cm",         1
    )

Then `AVAL <= 160` will lead to `AVALCAT1 == "<=160 cm"`, `AVAL`
in-between `160` and `170` will lead to `AVALCAT1 == "<=170 cm"`, and
`AVAL <= 170` will lead to `AVALCAT1 == ">170 cm"`.

However, we suggest to be more explicit when defining the `condition`,
to avoid overlap. In this case, the middle condition should be:
`AVAL <= 170 & AVAL > 160`

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_relative_flag.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)

## Examples

### Data setup

The following examples use the ADVS dataset below as a basis. It
contains vital signs data with some missing values (`NA`) that will
demonstrate how the function handles different scenarios.

    library(dplyr)
    library(tibble)

    advs <- tibble::tribble(
      ~USUBJID,       ~VSTEST,  ~AVAL,
      "01-701-1015", "Height", 147.32,
      "01-701-1015", "Weight",  53.98,
      "01-701-1023", "Height", 162.56,
      "01-701-1023", "Weight",     NA,
      "01-701-1028", "Height",     NA,
      "01-701-1028", "Weight",     NA,
      "01-701-1033", "Height", 175.26,
      "01-701-1033", "Weight",  88.45
    )

### Derive categorization variables without `by_vars`

In this example, we derive `AVALCAT1`, `AVALCA1N`, and `NEWCOL` without
using `by_vars`. The conditions must include all necessary filtering
logic, such as checking both `VSTEST` and `AVAL`. Records that don't
match any condition will have `NA` values for the new variables.

    definition <- exprs(
      ~condition,                        ~AVALCAT1, ~AVALCA1N,  ~NEWCOL,
      VSTEST == "Height" & AVAL > 160,   ">160 cm",         1, "extra1",
      VSTEST == "Height" & AVAL <= 160, "<=160 cm",         2, "extra2"
    )

    derive_vars_cat(
      dataset = advs,
      definition = definition
    )
    #> # A tibble: 8 × 6
    #>   USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N NEWCOL
    #>   <chr>       <chr>  <dbl> <chr>       <dbl> <chr>
    #> 1 01-701-1015 Height 147.  <=160 cm        2 extra2
    #> 2 01-701-1015 Weight  54.0 <NA>           NA <NA>
    #> 3 01-701-1023 Height 163.  >160 cm         1 extra1
    #> 4 01-701-1023 Weight  NA   <NA>           NA <NA>
    #> 5 01-701-1028 Height  NA   <NA>           NA <NA>
    #> 6 01-701-1028 Weight  NA   <NA>           NA <NA>
    #> 7 01-701-1033 Height 175.  >160 cm         1 extra1
    #> 8 01-701-1033 Weight  88.4 <NA>           NA <NA>  

### Derive categorization variables using `by_vars`

When using `by_vars`, the conditions are automatically scoped to records
matching each by group value. This simplifies the condition logic as you
don't need to include the by variable in each condition. Here we derive
categories for both Height and Weight measurements using
`by_vars = exprs(VSTEST)`.

    definition2 <- exprs(
      ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
      "Height",  AVAL > 160,  ">160 cm",         1,
      "Height", AVAL <= 160, "<=160 cm",         2,
      "Weight",   AVAL > 70,   ">70 kg",         1,
      "Weight",  AVAL <= 70,  "<=70 kg",         2
    )

    derive_vars_cat(
      dataset = advs,
      definition = definition2,
      by_vars = exprs(VSTEST)
    )
    #> # A tibble: 8 × 5
    #>   USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
    #>   <chr>       <chr>  <dbl> <chr>       <dbl>
    #> 1 01-701-1015 Height 147.  <=160 cm        2
    #> 2 01-701-1015 Weight  54.0 <=70 kg         2
    #> 3 01-701-1023 Height 163.  >160 cm         1
    #> 4 01-701-1023 Weight  NA   <NA>           NA
    #> 5 01-701-1028 Height  NA   <NA>           NA
    #> 6 01-701-1028 Weight  NA   <NA>           NA
    #> 7 01-701-1033 Height 175.  >160 cm         1
    #> 8 01-701-1033 Weight  88.4 >70 kg          1

### Using multiple conditions with explicit ranges

When you need more than two categories, you can define multiple
conditions. It's best practice to make conditions mutually exclusive
using explicit range definitions (e.g., `AVAL <= 170 & AVAL > 160`) to
avoid ambiguity, even though the function uses first-match logic.

    definition3 <- exprs(
      ~VSTEST,                ~condition,  ~AVALCAT1, ~AVALCA1N,
      "Height",               AVAL > 170,  ">170 cm",         1,
      "Height", AVAL <= 170 & AVAL > 160, "<=170 cm",         2,
      "Height",              AVAL <= 160, "<=160 cm",         3
    )

    derive_vars_cat(
      dataset = advs,
      definition = definition3,
      by_vars = exprs(VSTEST)
    )
    #> # A tibble: 8 × 5
    #>   USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
    #>   <chr>       <chr>  <dbl> <chr>       <dbl>
    #> 1 01-701-1015 Height 147.  <=160 cm        3
    #> 2 01-701-1015 Weight  54.0 <NA>           NA
    #> 3 01-701-1023 Height 163.  <=170 cm        2
    #> 4 01-701-1023 Weight  NA   <NA>           NA
    #> 5 01-701-1028 Height  NA   <NA>           NA
    #> 6 01-701-1028 Weight  NA   <NA>           NA
    #> 7 01-701-1033 Height 175.  >170 cm         1
    #> 8 01-701-1033 Weight  88.4 <NA>           NA

### Deriving categories based on reference ranges (`MCRITy` variables)

This example demonstrates deriving laboratory measurement criteria
variables (`MCRITyML` and `MCRITyMN`). The conditions use reference
range variables (like `ANRHI`) to create categories relative to normal
ranges, which is common in laboratory data analysis.

    adlb <- tibble::tribble(
      ~USUBJID,     ~PARAM, ~AVAL, ~AVALU,  ~ANRHI,
      "01-701-1015", "ALT",   150,  "U/L",      40,
      "01-701-1023", "ALT",    70,  "U/L",      40,
      "01-701-1036", "ALT",   130,  "U/L",      40,
      "01-701-1048", "ALT",    30,  "U/L",      40,
      "01-701-1015", "AST",    50,  "U/L",      35
    )

    definition_mcrit <- exprs(
      ~PARAM,                      ~condition,    ~MCRIT1ML, ~MCRIT1MN,
      "ALT",                    AVAL <= ANRHI,    "<=ANRHI",         1,
      "ALT", ANRHI < AVAL & AVAL <= 3 * ANRHI, ">1-3*ANRHI",         2,
      "ALT",                 3 * ANRHI < AVAL,   ">3*ANRHI",         3
    )

    adlb %>%
      derive_vars_cat(
        definition = definition_mcrit,
        by_vars = exprs(PARAM)
      )
    #> # A tibble: 5 × 7
    #>   USUBJID     PARAM  AVAL AVALU ANRHI MCRIT1ML   MCRIT1MN
    #>   <chr>       <chr> <dbl> <chr> <dbl> <chr>         <dbl>
    #> 1 01-701-1015 ALT     150 U/L      40 >3*ANRHI          3
    #> 2 01-701-1023 ALT      70 U/L      40 >1-3*ANRHI        2
    #> 3 01-701-1036 ALT     130 U/L      40 >3*ANRHI          3
    #> 4 01-701-1048 ALT      30 U/L      40 <=ANRHI           1
    #> 5 01-701-1015 AST      50 U/L      35 <NA>             NA

### Handling missing values and partial by groups

When using `by_vars`, records that don't match any by group in the
`definition` will have `NA` for all derived variables. In this example,
records with `VSTEST == "Weight"` will have `NA` values because only
"Height" conditions are defined. Additionally, records with missing
`AVAL` will result in `NA` for the categorization variables since
conditions cannot be evaluated.

    definition4 <- exprs(
      ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
      "Height",  AVAL > 160,  ">160 cm",         1,
      "Height", AVAL <= 160, "<=160 cm",         2
    )

    derive_vars_cat(
      dataset = advs,
      definition = definition4,
      by_vars = exprs(VSTEST)
    )
    #> # A tibble: 8 × 5
    #>   USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
    #>   <chr>       <chr>  <dbl> <chr>       <dbl>
    #> 1 01-701-1015 Height 147.  <=160 cm        2
    #> 2 01-701-1015 Weight  54.0 <NA>           NA
    #> 3 01-701-1023 Height 163.  >160 cm         1
    #> 4 01-701-1023 Weight  NA   <NA>           NA
    #> 5 01-701-1028 Height  NA   <NA>           NA
    #> 6 01-701-1028 Weight  NA   <NA>           NA
    #> 7 01-701-1033 Height 175.  >160 cm         1
    #> 8 01-701-1033 Weight  88.4 <NA>           NA

### Deriving multiple categorization variables simultaneously

You can derive any number of categorization variables in a single call.
This example creates three different categorization schemes (`AVALCAT1`,
`AVALCAT2`, and `AVALCAT3`) with their corresponding numeric flags, all
from the same set of conditions.

    definition5 <- exprs(
      ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N, ~AVALCAT2,      ~AVALCA2N, ~AVALCAT3,
      "Height",  AVAL > 160,  ">160 cm",         1,    "Tall",              1,   "Group A",
      "Height", AVAL <= 160, "<=160 cm",         2,   "Short",              2,   "Group B",
      "Weight",   AVAL > 70,   ">70 kg",         3,   "Heavy",              3,   "Group C",
      "Weight",  AVAL <= 70,  "<=70 kg",         4,   "Light",              4,   "Group D"
    )

    derive_vars_cat(
      dataset = advs,
      definition = definition5,
      by_vars = exprs(VSTEST)
    )
    #> # A tibble: 8 × 8
    #>   USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N AVALCAT2 AVALCA2N AVALCAT3
    #>   <chr>       <chr>  <dbl> <chr>       <dbl> <chr>       <dbl> <chr>
    #> 1 01-701-1015 Height 147.  <=160 cm        2 Short           2 Group B
    #> 2 01-701-1015 Weight  54.0 <=70 kg         4 Light           4 Group D
    #> 3 01-701-1023 Height 163.  >160 cm         1 Tall            1 Group A
    #> 4 01-701-1023 Weight  NA   <NA>           NA <NA>           NA <NA>
    #> 5 01-701-1028 Height  NA   <NA>           NA <NA>           NA <NA>
    #> 6 01-701-1028 Weight  NA   <NA>           NA <NA>           NA <NA>
    #> 7 01-701-1033 Height 175.  >160 cm         1 Tall            1 Group A
    #> 8 01-701-1033 Weight  88.4 >70 kg          3 Heavy           3 Group C 
