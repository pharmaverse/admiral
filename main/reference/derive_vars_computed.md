# Adds Variable(s) Computed from the Analysis Value of one or more Parameters

Adds Variable(s) computed from the analysis value of one or more
parameters. It is expected that the value of the new variable is defined
by an expression using the analysis values of other parameters, such as
addition/sum, subtraction/difference, multiplication/product,
division/ratio, exponentiation/logarithm, or by formula.  
  
For example Body Mass Index at Baseline (`BMIBL`) in `ADSL` can be
derived from of HEIGHT and WEIGHT parameters in `ADVS`.

## Usage

``` r
derive_vars_computed(
  dataset,
  dataset_add,
  by_vars,
  parameters,
  new_vars,
  filter_add = NULL,
  constant_by_vars = NULL,
  constant_parameters = NULL
)
```

## Arguments

- dataset:

  The variables specified by the `by_vars` parameter are expected.

  Default value

  :   none

- dataset_add:

  Additional dataset

  The variables specified by the `by_vars` parameter are expected.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the additional dataset after restricting it by the filter condition
  (`filter_add` parameter) and to the parameters specified by
  `parameters`.

  Default value

  :   none

- by_vars:

  Grouping variables

  Grouping variables uniquely identifying a set of records for which
  `new_vars` are to be calculated.

  Permitted values

  :   list of variables created by exprs()

  Default value

  :   none

- parameters:

  Required parameter codes

  It is expected that all parameter codes (`PARAMCD`) which are required
  to derive the new variable are specified for this parameter or the
  `constant_parameters` parameter.

  If observations should be considered which do not have a parameter
  code, e.g., if an SDTM dataset is used, temporary parameter codes can
  be derived by specifying a list of expressions. The name of the
  element defines the temporary parameter code and the expression
  defines the condition for selecting the records. For example,
  `parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the
  observations with `VSTESTCD == "HEIGHT"` from the input data
  (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
  observations, and adds them to the observations to consider.

  Unnamed elements in the list of expressions are considered as
  parameter codes. For example,
  `parameters = exprs(WEIGHT, HGHT = VSTESTCD == "HEIGHT")` uses the
  parameter code `"WEIGHT"` and creates a temporary parameter code
  `"HGHT"`.

  Permitted values

  :   A character vector of `PARAMCD` values or a list of expressions

  Default value

  :   none

- new_vars:

  Name of the newly created variables

  The specified variables are set to the specified values. The values of
  variables of the parameters specified by `parameters` can be accessed
  using `<variable name>.<parameter code>`. For example

      exprs(
        BMIBL = (AVAL.WEIGHT / (AVAL.HEIGHT/100)^2)
      )

  defines the value for the new variable.

  Variable names in the expression must not contain more than one dot.

  Permitted values

  :   List of variable-value pairs

  Default value

  :   none

- filter_add:

  Filter condition of additional dataset

  The specified condition is applied to the additional dataset before
  deriving the new variable, i.e., only observations fulfilling the
  condition are taken into account.

  Permitted values

  :   a condition

  Default value

  :   `NULL`

- constant_by_vars:

  By variables for constant parameters

  The constant parameters (parameters that are measured only once) are
  merged to the other parameters using the specified variables. (Refer
  to the Example)

  Permitted values

  :   list of variables

  Default value

  :   `NULL`

- constant_parameters:

  Required constant parameter codes

  It is expected that all the parameter codes (`PARAMCD`) which are
  required to derive the new variable and are measured only once are
  specified here. For example if BMI should be derived and height is
  measured only once while weight is measured at each visit. Height
  could be specified in the `constant_parameters` parameter. (Refer to
  the Example)

  If observations should be considered which do not have a parameter
  code, e.g., if an SDTM dataset is used, temporary parameter codes can
  be derived by specifying a list of expressions. The name of the
  element defines the temporary parameter code and the expression
  defines the condition for selecting the records. For example
  `constant_parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the
  observations with `VSTESTCD == "HEIGHT"` from the input data
  (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
  observations, and adds them to the observations to consider.

  Unnamed elements in the list of expressions are considered as
  parameter codes. For example,
  `constant_parameters = exprs(WEIGHT, HGHT = VSTESTCD == "HEIGHT")`
  uses the parameter code `"WEIGHT"` and creates a temporary parameter
  code `"HGHT"`.

  Permitted values

  :   A character vector of `PARAMCD` values or a list of expressions

  Default value

  :   `NULL`

## Value

The input dataset with the new variables added.

## Details

For each group (with respect to the variables specified for the
`by_vars` argument), the values of the new variables (`new_vars`) are
computed based on the parameters in the additional dataset
(`dataset_add`) and then the new variables are merged to the input
dataset (`dataset`).

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_cat.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)

## Examples

``` r
library(tibble)
library(dplyr)

# Example 1: Derive BMIBL
adsl <- tribble(
  ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01", "01-1302",   61,   "YEARS",
  "PILOT01", "17-1344",   64,   "YEARS"
)

advs <- tribble(
  ~STUDYID,  ~USUBJID,  ~PARAMCD, ~PARAM,        ~VISIT,      ~AVAL, ~AVALU, ~ABLFL,
  "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm",   "Y",
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE",   82.1, "kg",   "Y",
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2",    81.19, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4",    82.56, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6",    80.74, "kg",   NA,
  "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm",   "Y",
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE",  58.06, "kg",   "Y",
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2",    58.97, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4",    57.97, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6",    58.97, "kg",   NA
)

derive_vars_computed(
  dataset = adsl,
  dataset_add = advs,
  by_vars = exprs(STUDYID, USUBJID),
  parameters = c("WEIGHT", "HEIGHT"),
  new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
  filter_add = ABLFL == "Y"
)
#> # A tibble: 2 × 5
#>   STUDYID USUBJID   AGE AGEU  BMIBL
#>   <chr>   <chr>   <dbl> <chr> <dbl>
#> 1 PILOT01 01-1302    61 YEARS  26.0
#> 2 PILOT01 17-1344    64 YEARS  21.7
```
