# Add New Variable(s) to the Input Dataset Based on Variables from Another Dataset

Add new variable(s) to the input dataset based on variables from another
dataset. The observations to merge can be selected by a condition
(`filter_add` argument) and/or selecting the first or last observation
for each by group (`order` and `mode` argument).

## Usage

``` r
derive_vars_merged(
  dataset,
  dataset_add,
  by_vars,
  order = NULL,
  new_vars = NULL,
  filter_add = NULL,
  mode = NULL,
  exist_flag = NULL,
  true_value = "Y",
  false_value = NA_character_,
  missing_values = NULL,
  check_type = "warning",
  duplicate_msg = NULL,
  relationship = NULL
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

  The variables specified by the `by_vars`, the `new_vars`, and the
  `order` argument are expected.

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
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
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
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
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
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

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
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

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

- relationship:

  Expected merge-relationship between the `by_vars` variable(s) in
  `dataset` (input dataset) and the `dataset_add` (additional dataset)
  containing the additional `new_vars`.

  This argument is passed to the
  [`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  function. See
  <https://dplyr.tidyverse.org/reference/mutate-joins.html#arguments>
  for more details.

  Permitted values

  :   `"one-to-one"`, `"many-to-one"`

  Default value

  :   `NULL`

## Value

The output dataset contains all observations and variables of the input
dataset and additionally the variables specified for `new_vars` from the
additional dataset (`dataset_add`).

## Details

1.  The new variables (`new_vars`) are added to the additional dataset
    (`dataset_add`).

2.  The records from the additional dataset (`dataset_add`) are
    restricted to those matching the `filter_add` condition.

3.  If `order` is specified, for each by group the first or last
    observation (depending on `mode`) is selected.

4.  The variables specified for `new_vars` are merged to the input
    dataset using
    [`left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html).
    I.e., the output dataset contains all observations from the input
    dataset. For observations without a matching observation in the
    additional dataset the new variables are set as specified by
    `missing_values` (or to `NA` for variables not in `missing_values`).
    Observations in the additional dataset which have no matching
    observation in the input dataset are ignored.

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
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined_summary.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)

## Examples

### Note on usage versus [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)

The question between using `derive_vars_merged()` or the more powerful
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
comes down to how you need to select the observations to be merged.

- If the observations from `dataset_add` to merge can be selected by a
  condition (`filter_add`) using *only* variables from `dataset_add`,
  then always use `derive_vars_merged()` as it requires less resources
  (time and memory). A common example of this would be a randomization
  date in `ADSL`, where you are simply merging on a date from `DS`
  according to a certain `DSDECOD` condition such as
  `DSDECOD == "RANDOMIZATION"`.

- However, if the selection of the observations from `dataset_add` can
  depend on variables from *both* datasets, then use
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md).
  An example of this would be assigning period variables from `ADSL` to
  an `ADAE`, where you now need to check each adverse event start date
  against the period start and end dates to decide which period value to
  join.

### Basic merge of a full dataset

Merge all demographic variables onto a vital signs dataset.

- The variable `DOMAIN` exists in both datasets so note the use of
  `select(dm, -DOMAIN)` in the `dataset_add` argument. Without this an
  error would be issued to notify the user.

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)
    vs <- tribble(
      ~DOMAIN,  ~USUBJID, ~VSTESTCD, ~VISIT,      ~VSSTRESN, ~VSDTC,
      "VS",     "01",     "HEIGHT",  "SCREENING",     178.0, "2013-08-20",
      "VS",     "01",     "WEIGHT",  "SCREENING",      81.9, "2013-08-20",
      "VS",     "01",     "WEIGHT",  "BASELINE",       82.1, "2013-08-29",
      "VS",     "01",     "WEIGHT",  "WEEK 2",         81.9, "2013-09-15",
      "VS",     "01",     "WEIGHT",  "WEEK 4",         82.6, "2013-09-24",
      "VS",     "02",     "WEIGHT",  "BASELINE",       58.6, "2014-01-11"
    ) %>%
      mutate(STUDYID = "AB42")

    dm <- tribble(
      ~DOMAIN, ~USUBJID, ~AGE, ~AGEU,
      "DM",    "01",       61, "YEARS",
      "DM",    "02",       64, "YEARS",
      "DM",    "03",       85, "YEARS"
    ) %>%
      mutate(STUDYID = "AB42")

    derive_vars_merged(
      vs,
      dataset_add = select(dm, -DOMAIN),
      by_vars = exprs(STUDYID, USUBJID)
    ) %>%
      select(USUBJID, VSTESTCD, VISIT, VSSTRESN, AGE, AGEU)
    #> # A tibble: 6 × 6
    #>   USUBJID VSTESTCD VISIT     VSSTRESN   AGE AGEU
    #>   <chr>   <chr>    <chr>        <dbl> <dbl> <chr>
    #> 1 01      HEIGHT   SCREENING    178      61 YEARS
    #> 2 01      WEIGHT   SCREENING     81.9    61 YEARS
    #> 3 01      WEIGHT   BASELINE      82.1    61 YEARS
    #> 4 01      WEIGHT   WEEK 2        81.9    61 YEARS
    #> 5 01      WEIGHT   WEEK 4        82.6    61 YEARS
    #> 6 02      WEIGHT   BASELINE      58.6    64 YEARS

### Merge only the first/last value (`order` and `mode`)

Merge the last occurring weight for each subject to the demographics
dataset.

- To enable sorting by visit date
  [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md)
  is used to convert to a datetime, within the `order` argument.

- Then the `mode` argument is set to `"last"` to ensure the last sorted
  value is taken. Be cautious if `NA` values are possible in the `order`
  variables - see [Sort
  Order](https://pharmaverse.github.io/admiral/articles/generic.html#sort_order).

- The `filter_add` argument is used to restrict the vital signs records
  only to weight assessments.

    derive_vars_merged(
      dm,
      dataset_add = vs,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(convert_dtc_to_dtm(VSDTC)),
      mode = "last",
      new_vars = exprs(LSTWT = VSSTRESN),
      filter_add = VSTESTCD == "WEIGHT"
    ) %>%
      select(USUBJID, AGE, AGEU, LSTWT)
    #> # A tibble: 3 × 4
    #>   USUBJID   AGE AGEU  LSTWT
    #>   <chr>   <dbl> <chr> <dbl>
    #> 1 01         61 YEARS  82.6
    #> 2 02         64 YEARS  58.6
    #> 3 03         85 YEARS  NA  

### Handling duplicates (`check_type`)

The source records are checked regarding duplicates with respect to the
by variables and the order specified. By default, a warning is issued if
any duplicates are found. Note the results here with a new vital signs
dataset containing a duplicate last weight assessment date.

    vs_dup <- tribble(
      ~DOMAIN,  ~USUBJID, ~VSTESTCD, ~VISIT,      ~VSSTRESN, ~VSDTC,
      "VS",     "01",     "WEIGHT",  "WEEK 2",        81.1, "2013-09-24",
      "VS",     "01",     "WEIGHT",  "WEEK 4",        82.6, "2013-09-24"
    ) %>%
      mutate(STUDYID = "AB42")

    derive_vars_merged(
      dm,
      dataset_add = vs_dup,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(convert_dtc_to_dtm(VSDTC)),
      mode = "last",
      new_vars = exprs(LSTWT = VSSTRESN),
      filter_add = VSTESTCD == "WEIGHT"
    ) %>%
      select(USUBJID, AGE, AGEU, LSTWT)
    #> # A tibble: 3 × 4
    #>   USUBJID   AGE AGEU  LSTWT
    #>   <chr>   <dbl> <chr> <dbl>
    #> 1 01         61 YEARS  82.6
    #> 2 02         64 YEARS  NA
    #> 3 03         85 YEARS  NA
    #> Warning: Dataset contains duplicate records with respect to `STUDYID`, `USUBJID`, and
    #> `convert_dtc_to_dtm(VSDTC)`
    #> i Run `admiral::get_duplicates_dataset()` to access the duplicate records

For investigating the issue, the dataset of the duplicate source records
can be obtained by calling
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/get_duplicates_dataset.md):

    get_duplicates_dataset()
    #> Duplicate records with respect to `STUDYID`, `USUBJID`, and
    #> `convert_dtc_to_dtm(VSDTC)`.
    #> # A tibble: 2 × 9
    #>   STUDYID USUBJID convert_dtc_to_dtm(VSDT…¹ DOMAIN VSTESTCD VISIT VSSTRESN VSDTC
    #> * <chr>   <chr>   <dttm>                    <chr>  <chr>    <chr>    <dbl> <chr>
    #> 1 AB42    01      2013-09-24 00:00:00       VS     WEIGHT   WEEK…     81.1 2013…
    #> 2 AB42    01      2013-09-24 00:00:00       VS     WEIGHT   WEEK…     82.6 2013…
    #> # i abbreviated name: ¹`convert_dtc_to_dtm(VSDTC)`
    #> # i 1 more variable: LSTWT <dbl>

Common options to solve the issue:

- Specifying additional variables for `order` - this is the most common
  approach, adding something like a sequence variable.

- Restricting the source records by specifying/updating the `filter_add`
  argument.

- Setting `check_type = "none"` to ignore any duplicates, but then in
  this case the last occurring record would be chosen according to the
  sort order of the input `dataset_add`. This is not often advisable,
  unless the order has no impact on the result, as the temporary sort
  order can be prone to variation across an ADaM script.

### Modify values dependent on the merge (`new_vars` and `missing_values`)

For the last occurring weight for each subject, add a categorization of
which visit it occurred at to the demographics dataset.

- In the `new_vars` argument, other functions can be utilized to modify
  the merged values. For example, in the below case we want to
  categorize the visit as `"BASELINE"` or `"POST-BASELINE"` using
  [`if_else()`](https://dplyr.tidyverse.org/reference/if_else.html).

- The `missing_values` argument assigns a specific value for subjects
  with no matching observations - see subject `"03"` in the below
  example.

    derive_vars_merged(
      dm,
      dataset_add = vs,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(convert_dtc_to_dtm(VSDTC)),
      mode = "last",
      new_vars = exprs(
        LSTWTCAT = if_else(VISIT == "BASELINE", "BASELINE", "POST-BASELINE")
      ),
      filter_add = VSTESTCD == "WEIGHT",
      missing_values = exprs(LSTWTCAT = "MISSING")
    ) %>%
      select(USUBJID, AGE, AGEU, LSTWTCAT)
    #> # A tibble: 3 × 4
    #>   USUBJID   AGE AGEU  LSTWTCAT
    #>   <chr>   <dbl> <chr> <chr>
    #> 1 01         61 YEARS POST-BASELINE
    #> 2 02         64 YEARS BASELINE
    #> 3 03         85 YEARS MISSING      

### Check existence of records to merge (`exist_flag`, `true_value` and `false_value`)

Similar to the above example, now we prefer to have a separate flag
variable to show whether a selected record was merged.

- The name of the new variable is set with the `exist_flag` argument.

- The values of this new variable are assigned via the `true_value` and
  `false_value` arguments.

    derive_vars_merged(
      dm,
      dataset_add = vs,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(convert_dtc_to_dtm(VSDTC)),
      mode = "last",
      new_vars = exprs(
        LSTWTCAT = if_else(VISIT == "BASELINE", "BASELINE", "POST-BASELINE")
      ),
      filter_add = VSTESTCD == "WEIGHT",
      exist_flag = WTCHECK,
      true_value = "Y",
      false_value = "MISSING"
    ) %>%
      select(USUBJID, AGE, AGEU, LSTWTCAT, WTCHECK)
    #> # A tibble: 3 × 5
    #>   USUBJID   AGE AGEU  LSTWTCAT      WTCHECK
    #>   <chr>   <dbl> <chr> <chr>         <chr>
    #> 1 01         61 YEARS POST-BASELINE Y
    #> 2 02         64 YEARS BASELINE      Y
    #> 3 03         85 YEARS <NA>          MISSING

### Creating more than one variable from the merge (`new_vars`)

Derive treatment start datetime and associated imputation flags.

- In this example we first impute exposure datetime and associated flag
  variables as a separate first step to be used in the `order` argument.

- In the `new_vars` arguments, you can see how both datetime and the
  date and time imputation flags are all merged in one call.

    ex <- tribble(
      ~DOMAIN, ~USUBJID, ~EXSTDTC,
      "EX",    "01",     "2013-08-29",
      "EX",    "01",     "2013-09-16",
      "EX",    "02",     "2014-01-11",
      "EX",    "02",     "2014-01-25"
    ) %>%
      mutate(STUDYID = "AB42")

    ex_ext <- derive_vars_dtm(
      ex,
      dtc = EXSTDTC,
      new_vars_prefix = "EXST",
      highest_imputation = "M"
    )

    derive_vars_merged(
      dm,
      dataset_add = ex_ext,
      by_vars = exprs(STUDYID, USUBJID),
      new_vars = exprs(TRTSDTM = EXSTDTM, TRTSDTF = EXSTDTF, TRTSTMF = EXSTTMF),
      order = exprs(EXSTDTM),
      mode = "first"
    ) %>%
      select(USUBJID, TRTSDTM, TRTSDTF, TRTSTMF)
    #> # A tibble: 3 × 4
    #>   USUBJID TRTSDTM             TRTSDTF TRTSTMF
    #>   <chr>   <dttm>              <chr>   <chr>
    #> 1 01      2013-08-29 00:00:00 <NA>    H
    #> 2 02      2014-01-11 00:00:00 <NA>    H
    #> 3 03      NA                  <NA>    <NA>   

### Further examples

Further example usages of this function can be found in the
[`vignette("generic")`](https:/pharmaverse.github.io/admiral/main/articles/generic.md).
