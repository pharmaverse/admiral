# Add Variables from an Additional Dataset Based on Conditions from Both Datasets

The function adds variables from an additional dataset to the input
dataset. The selection of the observations from the additional dataset
can depend on variables from both datasets. For example, add the lowest
value (nadir) before the current observation.

## Usage

``` r
derive_vars_joined(
  dataset,
  dataset_add,
  by_vars = NULL,
  order = NULL,
  new_vars = NULL,
  tmp_obs_nr_var = NULL,
  join_vars = NULL,
  join_type,
  filter_add = NULL,
  first_cond_lower = NULL,
  first_cond_upper = NULL,
  filter_join = NULL,
  mode = NULL,
  exist_flag = NULL,
  true_value = "Y",
  false_value = NA_character_,
  missing_values = NULL,
  check_type = "warning"
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

  Additional dataset

  The variables specified by the `by_vars`, the `new_vars`, the
  `join_vars`, and the `order` argument are expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  Grouping variables

  The two datasets are joined by the specified variables.

  Variables can be renamed by naming the element, i.e.
  `by_vars = exprs(<name in input dataset> = <name in additional dataset>)`,
  similar to the `dplyr` joins.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- order:

  Sort order

  If the argument is set to a non-null value, for each observation of
  the input dataset the first or last observation from the joined
  dataset is selected with respect to the specified order. The specified
  variables are expected in the additional dataset (`dataset_add`). If a
  variable is available in both `dataset` and `dataset_add`, the one
  from `dataset_add` is used for the sorting.

  If an expression is named, e.g.,
  `exprs(EXSTDT = convert_dtc_to_dt(EXSTDTC), EXSEQ)`, a corresponding
  variable (`EXSTDT`) is added to the additional dataset and can be used
  in the filter conditions (`filter_add`, `filter_join`) and for
  `join_vars` and `new_vars`. The variable is not included in the output
  dataset.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- tmp_obs_nr_var:

  Temporary observation number

  The specified variable is added to the input dataset (`dataset`) and
  the additional dataset (`dataset_add`). It is set to the observation
  number with respect to `order`. For each by group (`by_vars`) the
  observation number starts with `1`. If there is more than one record
  for specific values for `by_vars` and `order`, all records get the
  same observation number. By default, a warning (see `check_type`) is
  issued in this case. The variable can be used in the conditions
  (`filter_join`, `first_cond_upper`, `first_cond_lower`). It can also
  be used to select consecutive observations or the last observation.

  The variable is not included in the output dataset. To include it
  specify it for `new_vars`.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- join_vars:

  Variables to use from additional dataset

  Any extra variables required from the additional dataset for
  `filter_join` should be specified for this argument. Variables
  specified for `new_vars` do not need to be repeated for `join_vars`.
  If a specified variable exists in both the input dataset and the
  additional dataset, the suffix ".join" is added to the variable from
  the additional dataset.

  If an expression is named, e.g.,
  `exprs(EXTDT = convert_dtc_to_dt(EXSTDTC))`, a corresponding variable
  is added to the additional dataset and can be used in the filter
  conditions (`filter_add`, `filter_join`) and for `new_vars`. The
  variable is not included in the output dataset.

  The variables are not included in the output dataset.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- join_type:

  Observations to keep after joining

  The argument determines which of the joined observations are kept with
  respect to the original observation. For example, if
  `join_type = "after"` is specified all observations after the original
  observations are kept.

  For example for confirmed response or BOR in the oncology setting or
  confirmed deterioration in questionnaires the confirmatory assessment
  must be after the assessment. Thus `join_type = "after"` could be
  used.

  Whereas, sometimes you might allow for confirmatory observations to
  occur prior to the observation. For example, to identify AEs occurring
  on or after seven days before a COVID AE. Thus `join_type = "all"`
  could be used.

  Permitted values

  :   `"before"`, `"after"`, `"all"`

  Default value

  :   none

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations from `dataset_add` fulfilling the specified
  condition are joined to the input dataset. If the argument is not
  specified, all observations are joined.

  Variables created by `order` or `new_vars` arguments can be used in
  the condition.

  The condition can include summary functions like
  [`all()`](https://rdrr.io/r/base/all.html) or
  [`any()`](https://rdrr.io/r/base/any.html). The additional dataset is
  grouped by the by variables (`by_vars`).

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- first_cond_lower:

  Condition for selecting range of data (before)

  If this argument is specified, the other observations are restricted
  from the last observation before the current observation where the
  specified condition is fulfilled up to the current observation. If the
  condition is not fulfilled for any of the other observations, no
  observations are considered.

  This argument should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only from a
  certain observation before the current observation up to the current
  observation. For an example, see the "Examples" section below.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- first_cond_upper:

  Condition for selecting range of data (after)

  If this argument is specified, the other observations are restricted
  up to the first observation where the specified condition is
  fulfilled. If the condition is not fulfilled for any of the other
  observations, no observations are considered.

  This argument should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only up to
  the confirmation assessment. For an example, see the "Examples"
  section below.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- filter_join:

  Filter for the joined dataset

  The specified condition is applied to the joined dataset. Therefore
  variables from both datasets `dataset` and `dataset_add` can be used.

  Variables created by `order` or `new_vars` arguments can be used in
  the condition.

  The condition can include summary functions like
  [`all()`](https://rdrr.io/r/base/all.html) or
  [`any()`](https://rdrr.io/r/base/any.html). The joined dataset is
  grouped by the original observations.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

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

- exist_flag:

  Exist flag

  If the argument is specified (e.g., `exist_flag = FLAG`), the
  specified variable (e.g., `FLAG`) is added to the input dataset. This
  variable will be the value provided in `true_value` for all selected
  records from `dataset_add` which are merged into the input dataset,
  and the value provided in `false_value` otherwise.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- true_value:

  True value

  The value for the specified variable `exist_flag`, applicable to the
  first or last observation (depending on the mode) of each by group.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"Y"`

- false_value:

  False value

  The value for the specified variable `exist_flag`, NOT applicable to
  the first or last observation (depending on the mode) of each by
  group.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- missing_values:

  Values for non-matching observations

  For observations of the input dataset (`dataset`) which do not have a
  matching observation in the additional dataset (`dataset_add`) the
  values of the specified variables are set to the specified value. Only
  variables specified for `new_vars` can be specified for
  `missing_values`.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   `NULL`

- check_type:

  Check uniqueness?

  If `"message"`, `"warning"` or `"error"` is specified, the specified
  message is issued if the observations of the (restricted) joined
  dataset are not unique with respect to the by variables and the order.

  This argument is ignored if `order` is not specified. In this case an
  error is issued independent of `check_type` if the restricted joined
  dataset contains more than one observation for any of the observations
  of the input dataset.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variables specified for `new_vars` from the
additional dataset (`dataset_add`).

## Details

1.  The variables specified by `order` are added to the additional
    dataset (`dataset_add`).

2.  The variables specified by `join_vars` are added to the additional
    dataset (`dataset_add`).

3.  The records from the additional dataset (`dataset_add`) are
    restricted to those matching the `filter_add` condition.

4.  The input dataset and the (restricted) additional dataset are left
    joined by the grouping variables (`by_vars`). If no grouping
    variables are specified, a full join is performed.

5.  If `first_cond_lower` is specified, for each observation of the
    input dataset the joined dataset is restricted to observations from
    the first observation where `first_cond_lower` is fulfilled (the
    observation fulfilling the condition is included) up to the
    observation of the input dataset. If for an observation of the input
    dataset the condition is not fulfilled, the observation is removed.

    If `first_cond_upper` is specified, for each observation of the
    input dataset the joined dataset is restricted to observations up to
    the first observation where `first_cond_upper` is fulfilled (the
    observation fulfilling the condition is included). If for an
    observation of the input dataset the condition is not fulfilled, the
    observation is removed.

    For an example, see the "Examples" section below.

6.  The joined dataset is restricted by the `filter_join` condition.

7.  If `order` is specified, for each observation of the input dataset
    the first or last observation (depending on `mode`) is selected.

8.  The variables specified for `new_vars` are created (if requested)
    and merged to the input dataset. I.e., the output dataset contains
    all observations from the input dataset. For observations without a
    matching observation in the joined dataset the new variables are set
    as specified by `missing_values` (or to `NA` for variables not in
    `missing_values`). Observations in the additional dataset which have
    no matching observation in the input dataset are ignored.

**Note:** This function creates temporary datasets which may be much
bigger than the input datasets. If this causes memory issues, please try
setting the admiral option `save_memory` to `TRUE` (see
[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/set_admiral_options.md)).
This reduces the memory consumption but increases the run-time.

## See also

[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/filter_joined.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_computed.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_transposed.md)

## Examples

### Note on usage versus [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)

The question between using
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
or the more powerful `derive_vars_joined()` comes down to how you need
to select the observations to be merged.

- If the observations from `dataset_add` to merge can be selected by a
  condition (`filter_add`) using *only* variables from `dataset_add`,
  then always use
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md)
  as it requires less resources (time and memory). A common example of
  this would be a randomization date in `ADSL`, where you are simply
  merging on a date from `DS` according to a certain `DSDECOD` condition
  such as `DSDECOD == "RANDOMIZATION"`.

- However, if the selection of the observations from `dataset_add` can
  depend on variables from *both* datasets, then use
  `derive_vars_joined()`. An example of this would be assigning period
  variables from `ADSL` to an `ADAE`, where you now need to check each
  adverse event start date against the period start and end dates to
  decide which period value to join.

### Basic join based on a generic time window (`filter_join`)

Derive a visit based on where the study day falls according to a
scheduled set of time windows.

- The `filter_join` argument here can check conditions using variables
  from both the `dataset` and `dataset_add`, so the study day is
  compared to the start and end of the time window.

- As no grouping variables are assigned using the `by_vars` argument, a
  full join is performed keeping all variables from `dataset_add`.

    library(tibble)
    library(lubridate)
    library(dplyr, warn.conflicts = FALSE)
    library(tidyr, warn.conflicts = FALSE)

    adbds <- tribble(
      ~USUBJID, ~ADY, ~AVAL,
      "1",       -33,    11,
      "1",        -7,    10,
      "1",         1,    12,
      "1",         8,    12,
      "1",        15,     9,
      "1",        20,    14,
      "1",        24,    12,
      "2",        -1,    13,
      "2",        13,     8
    ) %>%
      mutate(STUDYID = "AB42")

    windows <- tribble(
      ~AVISIT,    ~AWLO, ~AWHI,
      "BASELINE",   -30,     1,
      "WEEK 1",       2,     7,
      "WEEK 2",       8,    15,
      "WEEK 3",      16,    22,
      "WEEK 4",      23,    30
    )

    derive_vars_joined(
      adbds,
      dataset_add = windows,
      join_type = "all",
      filter_join = AWLO <= ADY & ADY <= AWHI
    ) %>%
      select(USUBJID, ADY, AWLO, AWHI, AVISIT)
    #> # A tibble: 9 × 5
    #>   USUBJID   ADY  AWLO  AWHI AVISIT
    #>   <chr>   <dbl> <dbl> <dbl> <chr>
    #> 1 1         -33    NA    NA <NA>
    #> 2 1          -7   -30     1 BASELINE
    #> 3 1           1   -30     1 BASELINE
    #> 4 1           8     8    15 WEEK 2
    #> 5 1          15     8    15 WEEK 2
    #> 6 1          20    16    22 WEEK 3
    #> 7 1          24    23    30 WEEK 4
    #> 8 2          -1   -30     1 BASELINE
    #> 9 2          13     8    15 WEEK 2  

### Join only the lowest/highest value occurring within a condition (`filter_join`, `order` and `mode`)

Derive the nadir value for each observation (i.e. the lowest value
occurring before) by subject.

- Note how `dataset` and `dataset_add` are the same here, so we are
  joining a dataset with itself. This enables us to compare records
  within the dataset to each other.

- Now we use `by_vars` as we only want to perform the join by subject.

- To find the lowest value we use the `order` and `mode` arguments.

- We subsequently need to check `ADY` to only check assessments
  occurring before. As this is not included in `by_vars` or `order`, we
  have to ensure it also gets joined by adding to `join_vars`. Then in
  `filter_join` note how `ADY.join < ADY` is used as the same variable
  exists in both datasets, so the version from `dataset_add` has `.join`
  added.

- According to the `AVAL` sort order used there could be duplicates
  (e.g. see subject `"1"` records at day 1 and 8), but given we only
  need to join `AVAL` itself here it doesn't actually matter to us which
  exact record is taken. So, in this example, we silence the uniqueness
  check by using `check_type = "none"`.

    derive_vars_joined(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(AVAL),
      new_vars = exprs(NADIR = AVAL),
      join_vars = exprs(ADY),
      join_type = "all",
      filter_join = ADY.join < ADY,
      mode = "first",
      check_type = "none"
    ) %>%
      select(USUBJID, ADY, AVAL, NADIR)
    #> # A tibble: 9 × 4
    #>   USUBJID   ADY  AVAL NADIR
    #>   <chr>   <dbl> <dbl> <dbl>
    #> 1 1         -33    11    NA
    #> 2 1          -7    10    11
    #> 3 1           1    12    10
    #> 4 1           8    12    10
    #> 5 1          15     9    10
    #> 6 1          20    14     9
    #> 7 1          24    12     9
    #> 8 2          -1    13    NA
    #> 9 2          13     8    13

### Filtering which records are joined from the additional dataset (`filter_add`)

Imagine we wanted to achieve the same as above, but we now want to
derive this allowing only post-baseline values to be possible for the
nadir.

- The `filter_add` argument can be used here as we only need to restrict
  the source data from `dataset_add`.

    derive_vars_joined(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(AVAL),
      new_vars = exprs(NADIR = AVAL),
      join_vars = exprs(ADY),
      join_type = "all",
      filter_add = ADY > 0,
      filter_join = ADY.join < ADY,
      mode = "first",
      check_type = "none"
    ) %>%
      select(USUBJID, ADY, AVAL, NADIR)
    #> # A tibble: 9 × 4
    #>   USUBJID   ADY  AVAL NADIR
    #>   <chr>   <dbl> <dbl> <dbl>
    #> 1 1         -33    11    NA
    #> 2 1          -7    10    NA
    #> 3 1           1    12    NA
    #> 4 1           8    12    12
    #> 5 1          15     9    12
    #> 6 1          20    14     9
    #> 7 1          24    12     9
    #> 8 2          -1    13    NA
    #> 9 2          13     8    NA

### Combining all of the above examples

Using all of the arguments demonstrated above, here is a more complex
example to add to `ADAE` the highest hemoglobin value occurring within
two weeks before each adverse event. Also join the day it occurred,
taking the earliest occurrence if more than one assessment with the same
value.

- Note how we used `mode = "last"` to get the highest lab value, but
  then as we wanted the earliest occurrence if more than one it means we
  need to add `desc(ADY)` to `order`. i.e. the last day when in
  descending order is the first.

    adae <- tribble(
      ~USUBJID, ~ASTDY,
      "1",           3,
      "1",          22,
      "2",           2
    ) %>%
      mutate(STUDYID = "AB42")

    adlb <- tribble(
      ~USUBJID, ~PARAMCD, ~ADY, ~AVAL,
      "1",      "HGB",       1,   8.5,
      "1",      "HGB",       3,   7.9,
      "1",      "HGB",       5,   8.9,
      "1",      "HGB",       8,   8.0,
      "1",      "HGB",       9,   8.0,
      "1",      "HGB",      16,   7.4,
      "1",      "ALB",       1,    42,
    ) %>%
      mutate(STUDYID = "AB42")

    derive_vars_joined(
      adae,
      dataset_add = adlb,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(AVAL, desc(ADY)),
      new_vars = exprs(HGB_MAX = AVAL, HGB_DY = ADY),
      join_type = "all",
      filter_add = PARAMCD == "HGB",
      filter_join = ASTDY - 14 <= ADY & ADY <= ASTDY,
      mode = "last"
    ) %>%
      select(USUBJID, ASTDY, HGB_MAX, HGB_DY)
    #> # A tibble: 3 × 4
    #>   USUBJID ASTDY HGB_MAX HGB_DY
    #>   <chr>   <dbl>   <dbl>  <dbl>
    #> 1 1           3     8.5      1
    #> 2 1          22     8        8
    #> 3 2           2    NA       NA

### Compute values in `new_vars` and `order`

Add to `ADAE` the number of days since the last dose of treatment, plus
1 day. If the dose occurs on the same day as the AE then include it as
the last dose.

- In the `new_vars` argument, other functions can be utilized to modify
  the joined values using variables from both `dataset` and
  `dataset_add`. For example, in the below case we want to calculate the
  number of days between the AE and the last dose using
  [`compute_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_duration.md).
  This function includes the plus 1 day as default.

- Also note how in this example `EXSDT` is created via the `order`
  argument and then used for `new_vars`, `filter_add` and `filter_join`.

- The reason to use `join_type = "all"` here instead of `"before"` is
  that we want to include any dose occurring on the same day as the AE,
  hence the `filter_join = EXSDT <= ASTDT`. Whereas using
  `join_type = "before"` would have resulted in the condition
  `EXSDT < ASTDT`. See the next example instead for
  `join_type = "before"`.

    adae <- tribble(
      ~USUBJID, ~ASTDT,
      "1",      "2020-02-02",
      "1",      "2020-02-04",
      "2",      "2021-01-08"
    ) %>%
      mutate(
        ASTDT = ymd(ASTDT),
        STUDYID = "AB42"
      )

    ex <- tribble(
      ~USUBJID, ~EXSDTC,
      "1",      "2020-01-10",
      "1",      "2020-01",
      "1",      "2020-01-20",
      "1",      "2020-02-03",
      "2",      "2021-01-05"
    ) %>%
      mutate(STUDYID = "AB42")

    derive_vars_joined(
      adae,
      dataset_add = ex,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(EXSDT = convert_dtc_to_dt(EXSDTC)),
      join_type = "all",
      new_vars = exprs(LDRELD = compute_duration(
        start_date = EXSDT, end_date = ASTDT
      )),
      filter_add = !is.na(EXSDT),
      filter_join = EXSDT <= ASTDT,
      mode = "last"
    ) %>%
      select(USUBJID, ASTDT, LDRELD)
    #> # A tibble: 3 × 3
    #>   USUBJID ASTDT      LDRELD
    #>   <chr>   <date>      <dbl>
    #> 1 1       2020-02-02     14
    #> 2 1       2020-02-04      2
    #> 3 2       2021-01-08      4

### Join records occurring before a condition (`join_type = "before"`)

In an arbitrary dataset where subjects have values of `"0"`, `"-"`,
`"+"` or `"++"`, for any value of `"0"` derive the last occurring `"++"`
day that occurs before the `"0"`.

- The `AVAL.join == "++"` in `filter_join`, along with `order` and
  `mode` taking the last day, identifies the target records to join from
  `dataset_add` for each observation of `dataset`.

- Then `join_type = "before"` is now used instead of
  `join_type = "all"`. This is because we only want to join the records
  occurring before the current observation in `dataset`. Including
  `AVAL == "0"` in `filter_join` ensures here that we only populate the
  new variable for records with `AVAL == "0"` in our `dataset`.

    myd <- tribble(
      ~USUBJID, ~ADY, ~AVAL,
      "1",         1, "++",
      "1",         2, "-",
      "1",         3, "0",
      "1",         4, "+",
      "1",         5, "++",
      "1",         6, "-",
      "2",         1, "-",
      "2",         2, "++",
      "2",         3, "+",
      "2",         4, "0",
      "2",         5, "-",
      "2",         6, "++",
      "2",         7, "0"
    ) %>%
      mutate(STUDYID = "AB42")

    derive_vars_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADY),
      mode = "last",
      new_vars = exprs(PREVPLDY = ADY),
      join_vars = exprs(AVAL),
      join_type = "before",
      filter_join = AVAL == "0" & AVAL.join == "++"
    ) %>%
      select(USUBJID, ADY, AVAL, PREVPLDY)
    #> # A tibble: 13 × 4
    #>    USUBJID   ADY AVAL  PREVPLDY
    #>    <chr>   <dbl> <chr>    <dbl>
    #>  1 1           1 ++          NA
    #>  2 1           2 -           NA
    #>  3 1           3 0            1
    #>  4 1           4 +           NA
    #>  5 1           5 ++          NA
    #>  6 1           6 -           NA
    #>  7 2           1 -           NA
    #>  8 2           2 ++          NA
    #>  9 2           3 +           NA
    #> 10 2           4 0            2
    #> 11 2           5 -           NA
    #> 12 2           6 ++          NA
    #> 13 2           7 0            6

### Join records occurring before a condition and checking all values in between (`first_cond_lower`, `join_type` and `filter_join`)

In the same example as above, now additionally check that in between the
`"++"` and the `"0"` all results must be either `"+"` or `"++"`.

- Firstly, `first_cond_lower = AVAL.join == "++"` is used so that for
  each observation of `dataset` the joined records from `dataset_add`
  are restricted to only include from the last occurring `"++"` before.
  This is necessary because of the use of a summary function in
  `filter_join` only on a subset of the joined observations as explained
  below.

- The `filter_join` condition used here now includes
  `all(AVAL.join %in% c("+", "++"))` to further restrict the joined
  records from `dataset_add` to only where all the values are either
  `"+"` or `"++"`.

- The `order` and `mode` arguments ensure only the day of the `"++"`
  value is joined. For example, for subject `"2"` it selects the day 2
  record instead of day 3, by using `"first"`.

    derive_vars_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADY),
      mode = "first",
      new_vars = exprs(PREVPLDY = ADY),
      join_vars = exprs(AVAL),
      join_type = "before",
      first_cond_lower = AVAL.join == "++",
      filter_join = AVAL == "0" & all(AVAL.join %in% c("+", "++"))
    ) %>%
      select(USUBJID, ADY, AVAL, PREVPLDY)
    #> # A tibble: 13 × 4
    #>    USUBJID   ADY AVAL  PREVPLDY
    #>    <chr>   <dbl> <chr>    <dbl>
    #>  1 1           1 ++          NA
    #>  2 1           2 -           NA
    #>  3 1           3 0           NA
    #>  4 1           4 +           NA
    #>  5 1           5 ++          NA
    #>  6 1           6 -           NA
    #>  7 2           1 -           NA
    #>  8 2           2 ++          NA
    #>  9 2           3 +           NA
    #> 10 2           4 0            2
    #> 11 2           5 -           NA
    #> 12 2           6 ++          NA
    #> 13 2           7 0            6

### Join records occurring after a condition checking all values in between (`first_cond_upper`, `join_type` and `filter_join`)

Similar to the above, now derive the first `"++"` day after any `"0"`
where all results in between are either `"+"` or `"++"`.

- Note how the main difference here is the use of `join_type = "after"`,
  `mode = "last"` and the `first_cond_upper` argument, instead of
  `first_cond_lower`.

    derive_vars_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADY),
      mode = "last",
      new_vars = exprs(NEXTPLDY = ADY),
      join_vars = exprs(AVAL),
      join_type = "after",
      first_cond_upper = AVAL.join == "++",
      filter_join = AVAL == "0" & all(AVAL.join %in% c("+", "++"))
    ) %>%
      select(USUBJID, ADY, AVAL, NEXTPLDY)
    #> # A tibble: 13 × 4
    #>    USUBJID   ADY AVAL  NEXTPLDY
    #>    <chr>   <dbl> <chr>    <dbl>
    #>  1 1           1 ++          NA
    #>  2 1           2 -           NA
    #>  3 1           3 0            5
    #>  4 1           4 +           NA
    #>  5 1           5 ++          NA
    #>  6 1           6 -           NA
    #>  7 2           1 -           NA
    #>  8 2           2 ++          NA
    #>  9 2           3 +           NA
    #> 10 2           4 0           NA
    #> 11 2           5 -           NA
    #> 12 2           6 ++          NA
    #> 13 2           7 0           NA

### Join a value from the next occurring record (`join_type = "after"`)

Add the value from the next occurring record as a new variable.

- The `join_type = "after"` here essentially acts as a lag to join
  variables from the next occurring record, and `mode = "first"` selects
  the first of these.

    derive_vars_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADY),
      mode = "first",
      new_vars = exprs(NEXTVAL = AVAL),
      join_vars = exprs(AVAL),
      join_type = "after"
    ) %>%
      select(USUBJID, ADY, AVAL, NEXTVAL)
    #> # A tibble: 13 × 4
    #>    USUBJID   ADY AVAL  NEXTVAL
    #>    <chr>   <dbl> <chr> <chr>
    #>  1 1           1 ++    -
    #>  2 1           2 -     0
    #>  3 1           3 0     +
    #>  4 1           4 +     ++
    #>  5 1           5 ++    -
    #>  6 1           6 -     <NA>
    #>  7 2           1 -     ++
    #>  8 2           2 ++    +
    #>  9 2           3 +     0
    #> 10 2           4 0     -
    #> 11 2           5 -     ++
    #> 12 2           6 ++    0
    #> 13 2           7 0     <NA>   

### Join records after a condition occurring in consecutive visits (`tmp_obs_nr_var`, `join_type` and `filter_join`)

Find the last occurring value on any of the next 3 unique visit days.

- The `tmp_obs_nr_var` argument can be useful as shown here to help pick
  out records happening before or after with respect to `order`, as you
  can see in the `filter_join`.

    derive_vars_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADY),
      mode = "last",
      new_vars = exprs(NEXTVAL = AVAL),
      tmp_obs_nr_var = tmp_obs_nr,
      join_vars = exprs(AVAL),
      join_type = "after",
      filter_join = tmp_obs_nr + 3 >= tmp_obs_nr.join
    ) %>%
      select(USUBJID, ADY, AVAL, NEXTVAL)
    #> # A tibble: 13 × 4
    #>    USUBJID   ADY AVAL  NEXTVAL
    #>    <chr>   <dbl> <chr> <chr>
    #>  1 1           1 ++    +
    #>  2 1           2 -     ++
    #>  3 1           3 0     -
    #>  4 1           4 +     -
    #>  5 1           5 ++    -
    #>  6 1           6 -     <NA>
    #>  7 2           1 -     0
    #>  8 2           2 ++    -
    #>  9 2           3 +     ++
    #> 10 2           4 0     0
    #> 11 2           5 -     0
    #> 12 2           6 ++    0
    #> 13 2           7 0     <NA>   

### Derive period variables (`APERIOD`, `APERSDT`, `APEREDT`)

Create a period reference dataset from `ADSL` and join this with `ADAE`
to identify within which period each AE occurred.

    adsl <- tribble(
      ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
      "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
      "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
    ) %>%
      mutate(across(ends_with("DT"), ymd)) %>%
      mutate(STUDYID = "AB42")

    period_ref <- create_period_dataset(
      adsl,
      new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)
    )

    period_ref
    #> # A tibble: 4 × 5
    #>   STUDYID USUBJID APERIOD APERSDT    APEREDT
    #>   <chr>   <chr>     <int> <date>     <date>
    #> 1 AB42    1             1 2021-01-04 2021-02-06
    #> 2 AB42    1             2 2021-02-07 2021-03-07
    #> 3 AB42    2             1 2021-02-02 2021-03-02
    #> 4 AB42    2             2 2021-03-03 2021-04-01

    adae <- tribble(
      ~USUBJID, ~ASTDT,
      "1",      "2021-01-01",
      "1",      "2021-01-05",
      "1",      "2021-02-05",
      "1",      "2021-03-05",
      "1",      "2021-04-05",
      "2",      "2021-02-15",
    ) %>%
      mutate(
        ASTDT = ymd(ASTDT),
        STUDYID = "AB42"
      )

    derive_vars_joined(
      adae,
      dataset_add = period_ref,
      by_vars = exprs(STUDYID, USUBJID),
      join_vars = exprs(APERSDT, APEREDT),
      join_type = "all",
      filter_join = APERSDT <= ASTDT & ASTDT <= APEREDT
    ) %>%
      select(USUBJID, ASTDT, APERSDT, APEREDT, APERIOD)
    #> # A tibble: 6 × 5
    #>   USUBJID ASTDT      APERSDT    APEREDT    APERIOD
    #>   <chr>   <date>     <date>     <date>       <int>
    #> 1 1       2021-01-01 NA         NA              NA
    #> 2 1       2021-01-05 2021-01-04 2021-02-06       1
    #> 3 1       2021-02-05 2021-01-04 2021-02-06       1
    #> 4 1       2021-03-05 2021-02-07 2021-03-07       2
    #> 5 1       2021-04-05 NA         NA              NA
    #> 6 2       2021-02-15 2021-02-02 2021-03-02       1

### Further examples

Further example usages of this function can be found in the
[`vignette("generic")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

Equivalent examples for using the `exist_flag`, `true_value`,
`false_value`, `missing_values` and `check_type` arguments can be found
in
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md).
