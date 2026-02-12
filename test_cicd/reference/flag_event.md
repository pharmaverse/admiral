# Create a `flag_event` Object

The `flag_event` object is used to define events as input for the
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md)
function.

## Usage

``` r
flag_event(dataset_name, condition = NULL, by_vars = NULL)
```

## Arguments

- dataset_name:

  Dataset name of the dataset to be used as input for the event. The
  name refers to the dataset specified for `source_datasets` in
  [`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md).

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- condition:

  Condition

  The condition is evaluated at the dataset referenced by
  `dataset_name`. For all by groups where it evaluates as `TRUE` at
  least once the new variable is set to the true value (`true_value`).

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- by_vars:

  Grouping variables

  If specified, the dataset is grouped by the specified variables before
  the condition is evaluated. If named elements are used in `by_vars`
  like `by_vars = exprs(USUBJID, EXLNKID = ECLNKID)`, the variables are
  renamed after the evaluation. If the `by_vars` element is not
  specified, the observations are grouped by the variables specified for
  the `by_vars` argument of
  [`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

## See also

[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_ef_msrc.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/test_cicd/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event_source.md),
[`query()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/tte_source.md)
