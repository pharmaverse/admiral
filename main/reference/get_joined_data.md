# Join Data for "joined" functions

The helper function joins the data for the "joined" functions. All
`.join` variables are included in the output dataset.

## Usage

``` r
get_joined_data(
  dataset,
  dataset_add,
  by_vars = NULL,
  join_vars = NULL,
  join_type,
  first_cond_lower = NULL,
  first_cond_upper = NULL,
  order = NULL,
  tmp_obs_nr_var = NULL,
  filter_add = NULL,
  filter_join = NULL,
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

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
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
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

- first_cond_lower:

  Condition for selecting range of data (before)

  If this argument is specified, the other observations are restricted
  from the first observation before the current observation where the
  specified condition is fulfilled up to the current observation. If the
  condition is not fulfilled for any of the other observations, no
  observations are considered, i.e., the observation is not flagged.

  This argument should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only from a
  certain observation before the current observation up to the current
  observation.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- first_cond_upper:

  Condition for selecting range of data (after)

  If this argument is specified, the other observations are restricted
  up to the first observation where the specified condition is
  fulfilled. If the condition is not fulfilled for any of the other
  observations, no observations are considered, i.e., the observation is
  not flagged.

  This argument should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only up to
  the confirmation assessment.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

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
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- tmp_obs_nr_var:

  Temporary observation number

  The specified variable is added to the input dataset (`dataset`) and
  the additional dataset (`dataset_add`). It is set to the observation
  number with respect to `order`. For each by group (`by_vars`) the
  observation number starts with `1`. The variable can be used in the
  conditions (`filter_join`, `first_cond_upper`, `first_cond_lower`). It
  can also be used to select consecutive observations or the last
  observation.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

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

5.  The joined dataset is restricted by the `filter_join` condition.
