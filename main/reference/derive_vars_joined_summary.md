# Summarize Variables from an Additional Dataset Based on Conditions from Both Datasets

The function summarizes variables from an additional dataset and adds
the summarized values as new variables to the input dataset. The
selection of the observations from the additional dataset can depend on
variables from both datasets. For example, all doses before the current
observation can be selected and the sum be added to the input dataset.

## Usage

``` r
derive_vars_joined_summary(
  dataset,
  dataset_add,
  by_vars = NULL,
  order = NULL,
  new_vars,
  tmp_obs_nr_var = NULL,
  join_vars = NULL,
  join_type,
  filter_add = NULL,
  first_cond_lower = NULL,
  first_cond_upper = NULL,
  filter_join = NULL,
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

  :   list of (optionally named) variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, ADY = ASTDY)`

  Default value

  :   `NULL`

- order:

  Sort order

  The specified variables are used to determine the order of the records
  if `first_cond_lower` or `first_cond_upper` is specified or if
  `join_type` equals `"before"` or `"after"`.

  If an expression is named, e.g.,
  `exprs(EXSTDT = convert_dtc_to_dt(EXSTDTC), EXSEQ)`, a corresponding
  variable (`EXSTDT`) is added to the additional dataset and can be used
  in the filter conditions (`filter_add`, `filter_join`) and for
  `join_vars` and `new_vars`. The variable is not included in the output
  dataset.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))` or `NULL`

  Default value

  :   `NULL`

- new_vars:

  Variables to add

  The new variables can be defined by named expressions, i.e.,
  `new_vars = exprs(<new variable> = <value>)`. The value must be
  defined such that it results in a single record per by group, e.g., by
  using a summary function like
  [`mean()`](https://rdrr.io/r/base/mean.html),
  [`sum()`](https://rdrr.io/r/base/sum.html), ...

  Permitted values

  :   list of named expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = "ml")`

  Default value

  :   none

- tmp_obs_nr_var:

  Temporary observation number

  The specified variable is added to the input dataset (`dataset`) and
  the restricted additional dataset (`dataset_add` after applying
  `filter_add`). It is set to the observation number with respect to
  `order`. For each by group (`by_vars`) the observation number starts
  with `1`. The variable can be used in the conditions (`filter_join`,
  `first_cond_upper`, `first_cond_lower`). It can also be used to select
  consecutive observations or the last observation.

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
  `exprs(EXSTDT = convert_dtc_to_dt(EXSTDTC))`, a corresponding variable
  is added to the additional dataset and can be used in the filter
  conditions (`filter_add`, `filter_join`) and for `new_vars`.

  The variables are not included in the output dataset.

  Permitted values

  :   list of variables or named expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(EXSTDY, EXSTDTM = convert_dtc_to_dtm(EXSTDTC))`

  Default value

  :   `NULL`

- join_type:

  Observations to keep after joining

  The argument determines which of the joined observations are kept with
  respect to the original observation. For example, if
  `join_type = "after"` is specified all observations after the original
  observations are kept.

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
  from the first observation before the current observation where the
  specified condition is fulfilled up to the current observation. If the
  condition is not fulfilled for any of the other observations, no
  observations are considered.

  This argument should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only from a
  certain observation before the current observation up to the current
  observation. For an example see the last example below.

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
  the confirmation assessment. For an example see the last example
  below.

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

- missing_values:

  Values for non-matching observations

  For observations of the input dataset (`dataset`) which do not have a
  matching observation in the additional dataset (`dataset_add`) the
  values of the specified variables are set to the specified value. Only
  variables specified for `new_vars` can be specified for
  `missing_values`.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   `NULL`

- check_type:

  Check uniqueness?

  If `"message"`, `"warning"` or `"error"` is specified, the specified
  message is issued if the observations of the input dataset (`dataset`)
  or the restricted additional dataset (`dataset_add` after applying
  `filter_add`) are not unique with respect to the by variables and the
  order.

  The uniqueness is checked only if `tmp_obs_nr_var`,
  `first_cond_lower`, or `first_cond_upper` is specified or `join_type`
  equals `"before"` or `"after"`.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variables specified for `new_vars` derived
from the additional dataset (`dataset_add`).

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

    For an example see the last example in the "Examples" section.

6.  The joined dataset is restricted by the `filter_join` condition.

7.  The variables specified for `new_vars` are created and merged to the
    input dataset. I.e., the output dataset contains all observations
    from the input dataset. For observations without a matching
    observation in the joined dataset the new variables are set as
    specified by `missing_values` (or to `NA` for variables not in
    `missing_values`). Observations in the additional dataset which have
    no matching observation in the input dataset are ignored.

**Note:** This function creates temporary datasets which may be much
bigger than the input datasets. If this causes memory issues, please try
setting the admiral option `save_memory` to `TRUE` (see
[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)).
This reduces the memory consumption but increases the run-time.

## See also

[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
[`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)

## Examples

The examples focus on the functionality specific to this function. For
examples of functionality common to all "joined" functions like
`filter_join`, `filter_add`, `join_vars`, ... please see the examples of
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md).

### Derive cumulative dose before event (`CUMDOSA`)

Deriving the cumulative actual dose up to the day of the adverse event
in the `ADAE` dataset.

- `USUBJID` is specified for `by_vars` to join the `ADAE` and the `ADEX`
  dataset by subject.

- `filter_join` is specified to restrict the `ADEX` dataset to the days
  up to the adverse event. `ADY.join` refers to the study day in `ADEX`.

- The new variable `CUMDOSA` is defined by the `new_vars` argument. It
  is set to the sum of `AVAL`.

- As `ADY` from `ADEX` is used in `filter_join` (but not in `new_vars`),
  it needs to be specified for `join_vars`.

- The `join_type` is set to `"all"` to consider all records in the
  joined dataset. `join_type = "before"` can't by used here because then
  doses at the same day as the adverse event would be excluded.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)

    adex <- tribble(
      ~USUBJID, ~ADY, ~AVAL,
      "1",         1,    10,
      "1",         8,    20,
      "1",        15,    10,
      "2",         8,     5
    )

    adae <- tribble(
      ~USUBJID, ~ADY, ~AEDECOD,
      "1",         2, "Fatigue",
      "1",         9, "Influenza",
      "1",        15, "Theft",
      "1",        15, "Fatigue",
      "2",         4, "Parasomnia",
      "3",         2, "Truancy"
    )

    derive_vars_joined_summary(
      dataset = adae,
      dataset_add = adex,
      by_vars = exprs(USUBJID),
      filter_join = ADY.join <= ADY,
      join_type = "all",
      join_vars = exprs(ADY),
      new_vars = exprs(CUMDOSA = sum(AVAL, na.rm = TRUE))
    )
    #> # A tibble: 6 × 4
    #>   USUBJID   ADY AEDECOD    CUMDOSA
    #>   <chr>   <dbl> <chr>        <dbl>
    #> 1 1           2 Fatigue         10
    #> 2 1           9 Influenza       30
    #> 3 1          15 Theft           40
    #> 4 1          15 Fatigue         40
    #> 5 2           4 Parasomnia      NA
    #> 6 3           2 Truancy         NA

### Define values for records without records in the additional dataset (`missing_values`)

By default, the new variables are set to `NA` for records without
matching records in the restricted additional dataset. This can be
changed by specifying the `missing_values` argument.

    derive_vars_joined_summary(
      dataset = adae,
      dataset_add = adex,
      by_vars = exprs(USUBJID),
      filter_join = ADY.join <= ADY,
      join_type = "all",
      join_vars = exprs(ADY),
      new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE)),
      missing_values = exprs(CUMDOSE = 0)
    )
    #> # A tibble: 6 × 4
    #>   USUBJID   ADY AEDECOD    CUMDOSE
    #>   <chr>   <dbl> <chr>        <dbl>
    #> 1 1           2 Fatigue         10
    #> 2 1           9 Influenza       30
    #> 3 1          15 Theft           40
    #> 4 1          15 Fatigue         40
    #> 5 2           4 Parasomnia       0
    #> 6 3           2 Truancy          0

### Selecting records (`join_type = "before"`, `join_type = "after"`)

The `join_type` argument can be used to select records from the
additional dataset. For example, if `join_type = "before"` is specified,
only records before the current observation are selected. If
`join_type = "after"` is specified, only records after the current
observation are selected.

To illustrate this, a variable (`SELECTED_DAYS`) is derived which
contains the selected days.

    mydata <- tribble(
      ~DAY,
      1,
      2,
      3,
      4,
      5
    )

    derive_vars_joined_summary(
      mydata,
      dataset_add = mydata,
      order = exprs(DAY),
      join_type = "before",
      new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
    )
    #> # A tibble: 5 × 2
    #>     DAY SELECTED_DAYS
    #>   <dbl> <chr>
    #> 1     1 <NA>
    #> 2     2 1
    #> 3     3 1, 2
    #> 4     4 1, 2, 3
    #> 5     5 1, 2, 3, 4

    derive_vars_joined_summary(
      mydata,
      dataset_add = mydata,
      order = exprs(DAY),
      join_type = "after",
      new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
    )
    #> # A tibble: 5 × 2
    #>     DAY SELECTED_DAYS
    #>   <dbl> <chr>
    #> 1     1 2, 3, 4, 5
    #> 2     2 3, 4, 5
    #> 3     3 4, 5
    #> 4     4 5
    #> 5     5 <NA>         

### Selecting records (`first_cond_lower`, `first_cond_upper`)

The `first_cond_lower` and `first_cond_upper` arguments can be used to
restrict the joined dataset to a certain range of records. For example,
if `first_cond_lower` is specified, the joined dataset is restricted to
the last observation before the current record where the condition is
fulfilled.

Please note:

- If the condition is not fulfilled for any of the records, no records
  are selected.

- The restriction implied by `join_type` is applied first.

- If a variable is contained in both `dataset` and `dataset_add` like
  `DAY` in the example below, `DAY` refers to the value from `dataset`
  and `DAY.join` to the value from `dataset_add`.

To illustrate this, a variable (`SELECTED_DAYS`) is derived which
contains the selected days.

    derive_vars_joined_summary(
      mydata,
      dataset_add = mydata,
      order = exprs(DAY),
      join_type = "before",
      first_cond_lower = DAY.join == 2,
      new_vars = exprs(SELECTED_DAYS = paste(sort(DAY), collapse = ", "))
    )
    #> # A tibble: 5 × 2
    #>     DAY SELECTED_DAYS
    #>   <dbl> <chr>
    #> 1     1 <NA>
    #> 2     2 <NA>
    #> 3     3 2
    #> 4     4 2, 3
    #> 5     5 2, 3, 4

    derive_vars_joined_summary(
      mydata,
      dataset_add = mydata,
      order = exprs(DAY),
      join_type = "after",
      first_cond_upper = DAY.join == 4,
      new_vars = exprs(SELECTED_DAYS = paste(DAY, collapse = ", "))
    )
    #> # A tibble: 5 × 2
    #>     DAY SELECTED_DAYS
    #>   <dbl> <chr>
    #> 1     1 2, 3, 4
    #> 2     2 3, 4
    #> 3     3 4
    #> 4     4 <NA>
    #> 5     5 <NA>

    derive_vars_joined_summary(
      mydata,
      dataset_add = mydata,
      order = exprs(DAY),
      join_type = "all",
      first_cond_lower = DAY.join == 2,
      first_cond_upper = DAY.join == 4,
      new_vars = exprs(SELECTED_DAYS = paste(sort(DAY), collapse = ", "))
    )
    #> # A tibble: 5 × 2
    #>     DAY SELECTED_DAYS
    #>   <dbl> <chr>
    #> 1     1 2, 3, 4
    #> 2     2 2, 3, 4
    #> 3     3 2, 3, 4
    #> 4     4 2, 3, 4
    #> 5     5 2, 3, 4      

### Derive weekly score if enough assessments are available

For each planned visit the average score within the week before the
visit should be derived if at least three assessments are available.

Please note that the condition for the number of assessments is
specified in `new_vars` and not in `filter_join`. This is because the
number of assessments within the week before the visit should be counted
but not the number of assessments available for the subject.

    planned_visits <- tribble(
      ~AVISIT,  ~ADY,
      "WEEK 1",    8,
      "WEEK 4",   29,
      "WEEK 8",   57
      ) %>%
      mutate(USUBJID = "1", .before = AVISIT)

    adqs <- tribble(
      ~ADY, ~AVAL,
         1,    10,
         2,    12,
         4,     9,
         5,     9,
         7,    10,
        25,    11,
        27,    10,
        29,    10,
        41,     8,
        42,     9,
        44,     5
    ) %>%
    mutate(USUBJID = "1")

    derive_vars_joined_summary(
      planned_visits,
      dataset_add = adqs,
      by_vars = exprs(USUBJID),
      filter_join = ADY - 7 <= ADY.join & ADY.join < ADY,
      join_type = "all",
      join_vars = exprs(ADY),
      new_vars = exprs(AVAL = if_else(n() >= 3, mean(AVAL, na.rm = TRUE), NA))
    )
    #> # A tibble: 3 × 4
    #>   USUBJID AVISIT   ADY  AVAL
    #>   <chr>   <chr>  <dbl> <dbl>
    #> 1 1       WEEK 1     8    10
    #> 2 1       WEEK 4    29    NA
    #> 3 1       WEEK 8    57    NA
