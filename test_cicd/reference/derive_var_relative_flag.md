# Flag Observations Before or After a Condition is Fulfilled

Flag all observations before or after the observation where a specified
condition is fulfilled for each by group. For example, the function
could be called to flag for each subject all observations before the
first disease progression or to flag all AEs after a specific AE.

## Usage

``` r
derive_var_relative_flag(
  dataset,
  by_vars,
  order,
  new_var,
  condition,
  mode,
  selection,
  inclusive,
  flag_no_ref_groups = TRUE,
  check_type = "warning"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Default value

  :   none

- by_vars:

  Grouping variables

  Default value

  :   none

- order:

  Sort order

  Within each by group the observations are ordered by the specified
  order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/test_cicd/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))`

  Default value

  :   none

- new_var:

  New variable

  The variable is added to the input dataset and set to `"Y"` for all
  observations before or after the condition is fulfilled. For all other
  observations it is set to `NA`.

  Default value

  :   none

- condition:

  Condition for Reference Observation

  The specified condition determines the reference observation. In the
  output dataset all observations before or after (`selection` argument)
  the reference observation are flagged.

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, for each by group the observations before
  or after (`selection` argument) the observation where the condition
  (`condition` argument) is fulfilled the *first* time is flagged in the
  output dataset. If `"last"` is specified, for each by group the
  observations before or after (`selection` argument) the observation
  where the condition (`condition` argument) is fulfilled the *last*
  time is flagged in the output dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- selection:

  Flag observations before or after the reference observation?

  Permitted values

  :   `"before"`, `"after"`

  Default value

  :   none

- inclusive:

  Flag the reference observation?

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   none

- flag_no_ref_groups:

  Should by groups without reference observation be flagged?

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the input dataset are not unique with
  respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

The input dataset with the new variable (`new_var`) added

## Details

For each by group (`by_vars` argument) the observations before or after
(`selection` argument) the observations where the condition (`condition`
argument) is fulfilled the first or last time (`order` argument and
`mode` argument) is flagged in the output dataset.

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_obs_number.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_transposed.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

# Flag all AEs after the first COVID AE
adae <- tribble(
  ~USUBJID, ~ASTDY, ~ACOVFL, ~AESEQ,
  "1",           2, NA,           1,
  "1",           5, "Y",          2,
  "1",           5, NA,           3,
  "1",          17, NA,           4,
  "1",          27, "Y",          5,
  "1",          32, NA,           6,
  "2",           8, NA,           1,
  "2",          11, NA,           2,
)

derive_var_relative_flag(
  adae,
  by_vars = exprs(USUBJID),
  order = exprs(ASTDY, AESEQ),
  new_var = PSTCOVFL,
  condition = ACOVFL == "Y",
  mode = "first",
  selection = "after",
  inclusive = FALSE,
  flag_no_ref_groups = FALSE
)
#> # A tibble: 8 × 5
#>   USUBJID ASTDY ACOVFL AESEQ PSTCOVFL
#>   <chr>   <dbl> <chr>  <dbl> <chr>   
#> 1 1           2 NA         1 NA      
#> 2 1           5 Y          2 NA      
#> 3 1           5 NA         3 Y       
#> 4 1          17 NA         4 Y       
#> 5 1          27 Y          5 Y       
#> 6 1          32 NA         6 Y       
#> 7 2           8 NA         1 NA      
#> 8 2          11 NA         2 NA      

response <- tribble(
  ~USUBJID, ~AVISITN, ~AVALC,
  "1",      0,        "PR",
  "1",      1,        "CR",
  "1",      2,        "CR",
  "1",      3,        "SD",
  "1",      4,        "NE",
  "2",      0,        "SD",
  "2",      1,        "PD",
  "2",      2,        "PD",
  "3",      0,        "SD",
  "4",      0,        "SD",
  "4",      1,        "PR",
  "4",      2,        "PD",
  "4",      3,        "SD",
  "4",      4,        "PR"
)

# Flag observations up to first PD for each patient
response %>%
  derive_var_relative_flag(
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    new_var = ANL02FL,
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )
#> # A tibble: 14 × 4
#>    USUBJID AVISITN AVALC ANL02FL
#>    <chr>     <dbl> <chr> <chr>  
#>  1 1             0 PR    Y      
#>  2 1             1 CR    Y      
#>  3 1             2 CR    Y      
#>  4 1             3 SD    Y      
#>  5 1             4 NE    Y      
#>  6 2             0 SD    Y      
#>  7 2             1 PD    Y      
#>  8 2             2 PD    NA     
#>  9 3             0 SD    Y      
#> 10 4             0 SD    Y      
#> 11 4             1 PR    Y      
#> 12 4             2 PD    Y      
#> 13 4             3 SD    NA     
#> 14 4             4 PR    NA     

# Flag observations up to first PD excluding baseline (AVISITN = 0) for each patient
response %>%
  restrict_derivation(
    derivation = derive_var_relative_flag,
    args = params(
      by_vars = exprs(USUBJID),
      order = exprs(AVISITN),
      new_var = ANL02FL,
      condition = AVALC == "PD",
      mode = "first",
      selection = "before",
      inclusive = TRUE
    ),
    filter = AVISITN > 0
  ) %>%
  arrange(USUBJID, AVISITN)
#> # A tibble: 14 × 4
#>    USUBJID AVISITN AVALC ANL02FL
#>    <chr>     <dbl> <chr> <chr>  
#>  1 1             0 PR    NA     
#>  2 1             1 CR    Y      
#>  3 1             2 CR    Y      
#>  4 1             3 SD    Y      
#>  5 1             4 NE    Y      
#>  6 2             0 SD    NA     
#>  7 2             1 PD    Y      
#>  8 2             2 PD    NA     
#>  9 3             0 SD    NA     
#> 10 4             0 SD    NA     
#> 11 4             1 PR    Y      
#> 12 4             2 PD    Y      
#> 13 4             3 SD    NA     
#> 14 4             4 PR    NA     
```
