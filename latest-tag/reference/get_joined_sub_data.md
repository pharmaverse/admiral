# Join Data for "joined" functions

The helper function joins the data for the "joined" functions. All
`.join` variables are included in the output dataset. It is called by
[`get_joined_data()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/get_joined_data.md)
to process each by group separately. This reduces the memory
consumption.

## Usage

``` r
get_joined_sub_data(
  dataset,
  dataset_add,
  by_vars,
  tmp_obs_nr_var,
  tmp_obs_nr_left,
  join_type,
  first_cond_upper,
  first_cond_lower,
  filter_join
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

- tmp_obs_nr_left:

  Temporary observation number for `dataset`

  The specified variable has to be in the input dataset (`dataset`) and
  has to be a unique key.

  Default value

  :   none

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

## Details

1.  The input dataset (`dataset`) and the additional dataset
    (`dataset_add`) are left joined by the grouping variables
    (`by_vars`). If no grouping variables are specified, a full join is
    performed.

2.  The joined dataset is restricted as specified by arguments
    `join_type`, `first_cond_upper`, and `first_cond_lower`. See
    argument descriptions for details.

3.  The joined dataset is restricted by the `filter_join` condition.
