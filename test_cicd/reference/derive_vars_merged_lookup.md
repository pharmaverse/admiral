# Merge Lookup Table with Source Dataset

Merge user-defined lookup table with the input dataset. Optionally print
a list of records from the input dataset that do not have corresponding
mapping from the lookup table.

## Usage

``` r
derive_vars_merged_lookup(
  dataset,
  dataset_add,
  by_vars,
  order = NULL,
  new_vars = NULL,
  mode = NULL,
  filter_add = NULL,
  check_type = "warning",
  duplicate_msg = NULL,
  print_not_mapped = TRUE
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_add:

  Lookup table

  The variables specified by the `by_vars` argument are expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  The input dataset and the selected observations from the additional
  dataset are merged by the specified variables.

  Variables can be renamed by naming the element, i.e.
  `by_vars = exprs(<name in input dataset> = <name in additional dataset>)`,
  similar to the `dplyr` joins.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- order:

  Sort order

  If the argument is set to a non-null value, for each by group the
  first or last observation from the additional dataset is selected with
  respect to the specified order.

  Variables defined by the `new_vars` argument can be used in the sort
  order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/test_cicd/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- new_vars:

  Variables to add

  The specified variables from the additional dataset are added to the
  output dataset. Variables can be renamed by naming the element, i.e.,
  `new_vars = exprs(<new name> = <old name>)`.

  For example `new_vars = exprs(var1, var2)` adds variables `var1` and
  `var2` from `dataset_add` to the input dataset.

  And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
  `old_var2` from `dataset_add` and adds them to the input dataset
  renaming `old_var2` to `new_var2`.

  Values of the added variables can be modified by specifying an
  expression. For example,
  `new_vars = LASTRSP = exprs(str_to_upper(AVALC))` adds the variable
  `LASTRSP` to the dataset and sets it to the upper case value of
  `AVALC`.

  If the argument is not specified or set to `NULL`, all variables from
  the additional dataset (`dataset_add`) are added. In the case when a
  variable exists in both datasets, an error is issued to ensure the
  user either adds to `by_vars`, removes or renames.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- mode:

  Selection mode

  Determines if the first or last observation is selected. If the
  `order` argument is specified, `mode` must be non-null.

  If the `order` argument is not specified, the `mode` argument is
  ignored.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   `NULL`

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations fulfilling the specified condition are taken into
  account for merging. If the argument is not specified, all
  observations are considered.

  Variables defined by the `new_vars` argument can be used in the filter
  condition.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- check_type:

  Check uniqueness?

  If `"warning"`, `"message"`, or `"error"` is specified, the specified
  message is issued if the observations of the (restricted) additional
  dataset are not unique with respect to the by variables and the order.

  If the `order` argument is not specified, the `check_type` argument is
  ignored: if the observations of the (restricted) additional dataset
  are not unique with respect to the by variables, an error is issued.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

- duplicate_msg:

  Message of unique check

  If the uniqueness check fails, the specified message is displayed.

  Permitted values

  :   a console message to be printed, e.g. `"Attention"` or for longer
      messages use `paste("Line 1", "Line 2")`

  Default value

  :   paste(
            "Dataset {.arg dataset_add} contains duplicate records with respect to",
            "{.var {vars2chr(by_vars)}}."
          )

- print_not_mapped:

  Print a list of unique `by_vars` values that do not have corresponding
  records from the lookup table?

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

## Value

The output dataset contains all observations and variables of the input
dataset, and add the variables specified in `new_vars` from the lookup
table specified in `dataset_add`. Optionally prints a list of unique
`by_vars` values that do not have corresponding records from the lookup
table (by specifying `print_not_mapped = TRUE`).

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_transposed.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)
vs <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID,        ~VISIT, ~VSTESTCD,       ~VSTEST,
  "PILOT01",    "VS", "01-1028",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "01-1028",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "01-1028",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "01-1028",      "WEEK 4",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "04-1325",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",      "WEEK 4",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "10-1027",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",      "WEEK 4",  "WEIGHT",      "Weight"
)

param_lookup <- tribble(
  ~VSTESTCD,                 ~VSTEST, ~PARAMCD,                       ~PARAM,
  "SYSBP", "Systolic Blood Pressure",  "SYSBP", "Syst Blood Pressure (mmHg)",
  "WEIGHT",                 "Weight", "WEIGHT",                "Weight (kg)",
  "HEIGHT",                 "Height", "HEIGHT",                "Height (cm)",
  "TEMP",              "Temperature",   "TEMP",            "Temperature (C)",
  "MAP",    "Mean Arterial Pressure",    "MAP",   "Mean Art Pressure (mmHg)",
  "BMI",           "Body Mass Index",    "BMI",    "Body Mass Index(kg/m^2)",
  "BSA",         "Body Surface Area",    "BSA",     "Body Surface Area(m^2)"
)

derive_vars_merged_lookup(
  dataset = vs,
  dataset_add = param_lookup,
  by_vars = exprs(VSTESTCD),
  new_vars = exprs(PARAMCD, PARAM),
  print_not_mapped = TRUE
)
#> All `VSTESTCD` are mapped.
#> # A tibble: 21 × 8
#>    STUDYID DOMAIN USUBJID VISIT       VSTESTCD VSTEST      PARAMCD PARAM        
#>    <chr>   <chr>  <chr>   <chr>       <chr>    <chr>       <chr>   <chr>        
#>  1 PILOT01 VS     01-1028 SCREENING   HEIGHT   Height      HEIGHT  Height (cm)  
#>  2 PILOT01 VS     01-1028 SCREENING   TEMP     Temperature TEMP    Temperature …
#>  3 PILOT01 VS     01-1028 BASELINE    TEMP     Temperature TEMP    Temperature …
#>  4 PILOT01 VS     01-1028 WEEK 4      TEMP     Temperature TEMP    Temperature …
#>  5 PILOT01 VS     01-1028 SCREENING 1 WEIGHT   Weight      WEIGHT  Weight (kg)  
#>  6 PILOT01 VS     01-1028 BASELINE    WEIGHT   Weight      WEIGHT  Weight (kg)  
#>  7 PILOT01 VS     01-1028 WEEK 4      WEIGHT   Weight      WEIGHT  Weight (kg)  
#>  8 PILOT01 VS     04-1325 SCREENING   HEIGHT   Height      HEIGHT  Height (cm)  
#>  9 PILOT01 VS     04-1325 SCREENING   TEMP     Temperature TEMP    Temperature …
#> 10 PILOT01 VS     04-1325 BASELINE    TEMP     Temperature TEMP    Temperature …
#> # ℹ 11 more rows
```
