# Filter Observations Taking Other Observations into Account

The function filters observation using a condition taking other
observations into account. For example, it could select all observations
with `AVALC == "Y"` and `AVALC == "Y"` for at least one subsequent
observation. The input dataset is joined with itself to enable
conditions taking variables from both the current observation and the
other observations into account. The suffix ".join" is added to the
variables from the subsequent observations.

An example usage might be checking if a patient received two required
medications within a certain timeframe of each other.

In the oncology setting, for example, we use such processing to check if
a response value can be confirmed by a subsequent assessment. This is
commonly used in endpoints such as best overall response.

## Usage

``` r
filter_joined(
  dataset,
  dataset_add,
  by_vars,
  join_vars,
  join_type,
  first_cond_lower = NULL,
  first_cond_upper = NULL,
  order = NULL,
  tmp_obs_nr_var = NULL,
  filter_add = NULL,
  filter_join,
  check_type = "warning"
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` and `order` arguments are
  expected to be in the dataset.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- dataset_add:

  Additional dataset

  The variables specified for `by_vars`, `join_vars`, and `order` are
  expected.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   none

- by_vars:

  By variables

  The specified variables are used as by variables for joining the input
  dataset with itself.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- join_vars:

  Variables to keep from joined dataset

  The variables needed from the other observations should be specified
  for this parameter. The specified variables are added to the joined
  dataset with suffix ".join". For example to select all observations
  with `AVALC == "Y"` and `AVALC == "Y"` for at least one subsequent
  visit `join_vars = exprs(AVALC, AVISITN)` and
  `filter_join = AVALC == "Y" & AVALC.join == "Y" & AVISITN < AVISITN.join`
  could be specified.

  The `*.join` variables are not included in the output dataset.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

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

- first_cond_lower:

  Condition for selecting range of data (before)

  If this argument is specified, the other observations are restricted
  from the first observation before the current observation where the
  specified condition is fulfilled up to the current observation. If the
  condition is not fulfilled for any of the other observations, no
  observations are considered, i.e., the observation is not flagged.

  This parameter should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only from a
  certain observation before the current observation up to the current
  observation. For examples see the "Examples" section below.

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

  This parameter should be specified if `filter_join` contains summary
  functions which should not apply to all observations but only up to
  the confirmation assessment. For examples see the "Examples" section
  below.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- order:

  Order

  The observations are ordered by the specified order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of expressions created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(ADT, desc(AVAL))` or `NULL`

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
  (`filter_join`, `first_cond_upper`, `first_cond_lower`). It is not
  included in the output dataset. It can also be used to select
  consecutive observations or the last observation (see example below).

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   `NULL`

- filter_add:

  Filter for additional dataset (`dataset_add`)

  Only observations from `dataset_add` fulfilling the specified
  condition are joined to the input dataset. If the argument is not
  specified, all observations are joined.

  Variables created by the `order` argument can be used in the
  condition.

  The condition can include summary functions. The additional dataset is
  grouped by the by variables (`by_vars`).

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- filter_join:

  Condition for selecting observations

  The filter is applied to the joined dataset for selecting the
  confirmed observations. The condition can include summary functions
  like [`all()`](https://rdrr.io/r/base/all.html) or
  [`any()`](https://rdrr.io/r/base/any.html). The joined dataset is
  grouped by the original observations. I.e., the summary function are
  applied to all observations up to the confirmation observation. For
  example in the oncology setting when using this function for confirmed
  best overall response,
  `filter_join = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) & count_vals(var = AVALC.join, val = "NE") <= 1`
  selects observations with response "CR" and for all observations up to
  the confirmation observation the response is "CR" or "NE" and there is
  at most one "NE".

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   none

- check_type:

  Check uniqueness?

  If `"message"`, `"warning"`, or `"error"` is specified, the specified
  message is issued if the observations of the input dataset are not
  unique with respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

A subset of the observations of the input dataset. All variables of the
input dataset are included in the output dataset.

## Details

The following steps are performed to produce the output dataset.

### Step 1

- The variables specified by `order` are added to the additional dataset
  (`dataset_add`).

- The variables specified by `join_vars` are added to the additional
  dataset (`dataset_add`).

- The records from the additional dataset (`dataset_add`) are restricted
  to those matching the `filter_add` condition.

Then the input dataset (`dataset`) is joined with the restricted
additional dataset by the variables specified for `by_vars`. From the
additional dataset only the variables specified for `join_vars` are
kept. The suffix ".join" is added to those variables which are also
present in the input dataset.

For example, for `by_vars = USUBJID`,
`join_vars = exprs(AVISITN, AVALC)` and input dataset and additional
dataset

    # A tibble: 2 x 4
    USUBJID AVISITN AVALC  AVAL
    <chr>     <dbl> <chr> <dbl>
    1             1 Y         1
    1             2 N         0

the joined dataset is

    A tibble: 4 x 6
    USUBJID AVISITN AVALC  AVAL AVISITN.join AVALC.join
    <chr>     <dbl> <chr> <dbl>        <dbl> <chr>
    1             1 Y         1            1 Y
    1             1 Y         1            2 N
    1             2 N         0            1 Y
    1             2 N         0            2 N

### Step 2

The joined dataset is restricted to observations with respect to
`join_type` and `order`.

The dataset from the example in the previous step with
`join_type = "after"` and `order = exprs(AVISITN)` is restricted to

    A tibble: 4 x 6
    USUBJID AVISITN AVALC  AVAL AVISITN.join AVALC.join
    <chr>     <dbl> <chr> <dbl>        <dbl> <chr>
    1             1 Y         1            2 N

### Step 3

If `first_cond_lower` is specified, for each observation of the input
dataset the joined dataset is restricted to observations from the first
observation where `first_cond_lower` is fulfilled (the observation
fulfilling the condition is included) up to the observation of the input
dataset. If for an observation of the input dataset the condition is not
fulfilled, the observation is removed.

If `first_cond_upper` is specified, for each observation of the input
dataset the joined dataset is restricted to observations up to the first
observation where `first_cond_upper` is fulfilled (the observation
fulfilling the condition is included). If for an observation of the
input dataset the condition is not fulfilled, the observation is
removed.

For an example see the last example in the "Examples" section.

### Step 4

The joined dataset is grouped by the observations from the input dataset
and restricted to the observations fulfilling the condition specified by
`filter_join`.

### Step 5

The first observation of each group is selected and the `*.join`
variables are dropped.

**Note:** This function creates temporary datasets which may be much
bigger than the input datasets. If this causes memory issues, please try
setting the admiral option `save_memory` to `TRUE` (see
[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)).
This reduces the memory consumption but increases the run-time.

## See also

[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md)

Utilities for Filtering Observations:
[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md),
[`filter_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_exist.md),
[`filter_extreme()`](https:/pharmaverse.github.io/admiral/main/reference/filter_extreme.md),
[`filter_not_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_not_exist.md),
[`filter_relative()`](https:/pharmaverse.github.io/admiral/main/reference/filter_relative.md),
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md),
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)

## Examples

### Filter records considering other records (`filter_join`, `join_vars`)

In this example, the input dataset should be restricted to records with
a duration longer than 30 and where a COVID AE (`ACOVFL == "Y"`)
occurred before or up to seven days after the record. The condition for
restricting the records is specified by the `filter_join` argument.
Variables from the other records are referenced by variable names with
the suffix `.join`. These variables have to be specified for the
`join_vars` argument. As records before *and* after the current record
should be considered, `join_type = "all"` is specified.

    library(tibble)

    adae <- tribble(
      ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
      "1",        10, "N",          1,
      "1",        21, "N",         50,
      "1",        23, "Y",         14,
      "1",        32, "N",         31,
      "1",        42, "N",         20,
      "2",        11, "Y",         13,
      "2",        23, "N",          2,
      "3",        13, "Y",         12,
      "4",        14, "N",         32,
      "4",        21, "N",         41
    )

    filter_joined(
      adae,
      dataset_add = adae,
      by_vars = exprs(USUBJID),
      join_vars = exprs(ACOVFL, ADY),
      join_type = "all",
      filter_join = ADURN > 30 & ACOVFL.join == "Y" & ADY.join <= ADY + 7
    )
    #> # A tibble: 2 × 4
    #>   USUBJID   ADY ACOVFL ADURN
    #>   <chr>   <dbl> <chr>  <dbl>
    #> 1 1          21 N         50
    #> 2 1          32 N         31

### Considering only records after the current one (`join_type = "after"`)

In this example, the input dataset is restricted to records with
`AVALC == "Y"` and `AVALC == "Y"` at a subsequent visit.
`join_type = "after"` is specified to consider only records after the
current one. Please note that the `order` argument must be specified, as
otherwise it is not possible to determine which records are after the
current record.

    data <- tribble(
      ~USUBJID, ~AVISITN, ~AVALC,
      "1",      1,        "Y",
      "1",      2,        "N",
      "1",      3,        "Y",
      "1",      4,        "N",
      "2",      1,        "Y",
      "2",      2,        "N",
      "3",      1,        "Y",
      "4",      1,        "N",
      "4",      2,        "N",
    )

    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC, AVISITN),
      join_type = "after",
      order = exprs(AVISITN),
      filter_join = AVALC == "Y" & AVALC.join == "Y"
    )
    #> # A tibble: 1 × 3
    #>   USUBJID AVISITN AVALC
    #>   <chr>     <dbl> <chr>
    #> 1 1             1 Y    

### Considering a range of records only (`first_cond_lower`, `first_cond_upper`)

Consider the following data.

    myd <- tribble(
      ~subj, ~day, ~val,
      "1",      1, "++",
      "1",      2, "-",
      "1",      3, "0",
      "1",      4, "+",
      "1",      5, "++",
      "1",      6, "-",
      "2",      1, "-",
      "2",      2, "++",
      "2",      3, "+",
      "2",      4, "0",
      "2",      5, "-",
      "2",      6, "++"
    )

To select `"0"` where all results from the first `"++"` before the `"0"`
up to the `"0"` (excluding the `"0"`) are `"+"` or `"++"` the
`first_cond_lower` argument and `join_type = "before"` are specified.

    filter_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(subj),
      order = exprs(day),
      join_vars = exprs(val),
      join_type = "before",
      first_cond_lower = val.join == "++",
      filter_join = val == "0" & all(val.join %in% c("+", "++"))
    )
    #> # A tibble: 1 × 3
    #>   subj    day val
    #>   <chr> <dbl> <chr>
    #> 1 2         4 0    

To select `"0"` where all results from the `"0"` (excluding the `"0"`)
up to the first `"++"` after the `"0"` are `"+"` or `"++"` the
`first_cond_upper` argument and `join_type = "after"` are specified.

    filter_joined(
      myd,
      dataset_add = myd,
      by_vars = exprs(subj),
      order = exprs(day),
      join_vars = exprs(val),
      join_type = "after",
      first_cond_upper = val.join == "++",
      filter_join = val == "0" & all(val.join %in% c("+", "++"))
    )
    #> # A tibble: 1 × 3
    #>   subj    day val
    #>   <chr> <dbl> <chr>
    #> 1 1         3 0    

### Considering only records up to a condition (`first_cond_upper`)

In this example from deriving confirmed response in oncology, the
records with

- `AVALC == "CR"`,

- `AVALC == "CR"` at a subsequent visit,

- only `"CR"` or `"NE"` in between, and

- at most one `"NE"` in between

should be selected. The other records to be considered are restricted to
those up to the first occurrence of `"CR"` by specifying the
`first_cond_upper` argument. The
[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md)
function is used to count the `"NE"`s for the last condition.

    data <- tribble(
      ~USUBJID, ~AVISITN, ~AVALC,
      "1",      1,        "PR",
      "1",      2,        "CR",
      "1",      3,        "NE",
      "1",      4,        "CR",
      "1",      5,        "NE",
      "2",      1,        "CR",
      "2",      2,        "PR",
      "2",      3,        "CR",
      "3",      1,        "CR",
      "4",      1,        "CR",
      "4",      2,        "NE",
      "4",      3,        "NE",
      "4",      4,        "CR",
      "4",      5,        "PR"
    )

    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      order = exprs(AVISITN),
      first_cond_upper = AVALC.join == "CR",
      filter_join = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) &
        count_vals(var = AVALC.join, val = "NE") <= 1
    )
    #> # A tibble: 1 × 3
    #>   USUBJID AVISITN AVALC
    #>   <chr>     <dbl> <chr>
    #> 1 1             2 CR   

### Considering order of values ([`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md), [`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md))

In this example from deriving confirmed response in oncology, records
with

- `AVALC == "PR"`,

- `AVALC == "CR"` or `AVALC == "PR"` at a subsequent visit at least 20
  days later,

- only `"CR"`, `"PR"`, or `"NE"` in between,

- at most one `"NE"` in between, and

- `"CR"` is not followed by `"PR"`

should be selected. The last condition is realized by using
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)
and
[`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md),
ensuring that the first occurrence of `"CR"` is after the last
occurrence of `"PR"`. The second call to
[`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md)
in the condition is required to cover the case of no `"CR"`s (the
[`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md)
call returns `NA` then).

    data <- tribble(
      ~USUBJID, ~ADY, ~AVALC,
      "1",         6, "PR",
      "1",        12, "CR",
      "1",        24, "NE",
      "1",        32, "CR",
      "1",        48, "PR",
      "2",         3, "PR",
      "2",        21, "CR",
      "2",        33, "PR",
      "3",        11, "PR",
      "4",         7, "PR",
      "4",        12, "NE",
      "4",        24, "NE",
      "4",        32, "PR",
      "4",        55, "PR"
    )

    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC, ADY),
      join_type = "after",
      order = exprs(ADY),
      first_cond_upper = AVALC.join %in% c("CR", "PR") & ADY.join - ADY >= 20,
      filter_join = AVALC == "PR" &
        all(AVALC.join %in% c("CR", "PR", "NE")) &
        count_vals(var = AVALC.join, val = "NE") <= 1 &
        (
          min_cond(var = ADY.join, cond = AVALC.join == "CR") >
            max_cond(var = ADY.join, cond = AVALC.join == "PR") |
            count_vals(var = AVALC.join, val = "CR") == 0
        )
    )
    #> # A tibble: 1 × 3
    #>   USUBJID   ADY AVALC
    #>   <chr>   <dbl> <chr>
    #> 1 4          32 PR   

### Considering the order of records (`tmp_obs_nr_var`)

In this example, the records with `CRIT1FL == "Y"` at two consecutive
visits or at the last visit should be selected. A temporary order
variable is created by specifying the `tmp_obs_nr_var` argument. Then it
is used in `filter_join`. The temporary variable doesn't need to be
specified for `join_vars`.

    data <- tribble(
      ~USUBJID, ~AVISITN, ~CRIT1FL,
      "1",      1,        "Y",
      "1",      2,        "N",
      "1",      3,        "Y",
      "1",      5,        "N",
      "2",      1,        "Y",
      "2",      3,        "Y",
      "2",      5,        "N",
      "3",      1,        "Y",
      "4",      1,        "Y",
      "4",      2,        "N",
    )

    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      tmp_obs_nr_var = tmp_obs_nr,
      join_vars = exprs(CRIT1FL),
      join_type = "all",
      order = exprs(AVISITN),
      filter_join = CRIT1FL == "Y" & CRIT1FL.join == "Y" &
        (tmp_obs_nr + 1 == tmp_obs_nr.join | tmp_obs_nr == max(tmp_obs_nr.join))
    )
    #> # A tibble: 2 × 3
    #>   USUBJID AVISITN CRIT1FL
    #>   <chr>     <dbl> <chr>
    #> 1 2             1 Y
    #> 2 3             1 Y      
