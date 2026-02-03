# Add a Variable Flagging the First or Last Observation Within Each By Group

Add a variable flagging the first or last observation within each by
group

## Usage

``` r
derive_var_extreme_flag(
  dataset,
  by_vars,
  order,
  new_var,
  mode,
  true_value = "Y",
  false_value = NA_character_,
  flag_all = FALSE,
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

- by_vars:

  Grouping variables

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- order:

  Sort order

  The first or last observation is determined with respect to the
  specified order.

  For handling of `NA`s in sorting variables see the "Sort Order"
  section in
  [`vignette("generic")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/generic.md).

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- new_var:

  Variable to add

  The specified variable is added to the output dataset. It is set to
  the value set in `true_value` for the first or last observation
  (depending on the mode) of each by group.

  Permitted values

  :   an unquoted symbol, e.g., `AVAL`

  Default value

  :   none

- mode:

  Flag mode

  Determines of the first or last observation is flagged.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- true_value:

  True value

  The value for the specified variable `new_var`, applicable to the
  first or last observation (depending on the mode) of each by group.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"Y"`

- false_value:

  False value

  The value for the specified variable `new_var`, NOT applicable to the
  first or last observation (depending on the mode) of each by group.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `NA_character_`

- flag_all:

  Flag setting

  A logical value where if set to `TRUE`, all records are flagged and no
  error or warning is issued if the first or last record is not unique.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- check_type:

  Check uniqueness?

  If `"warning"` or `"error"` is specified, the specified message is
  issued if the observations of the input dataset are not unique with
  respect to the by variables and the order.

  Permitted values

  :   `"none"`, `"message"`, `"warning"`, `"error"`

  Default value

  :   `"warning"`

## Value

The input dataset with the new flag variable added

## Details

For each group (with respect to the variables specified for the
`by_vars` parameter), `new_var` is set to `"Y"` for the first or last
observation (with respect to the order specified for the `order`
parameter and the flag mode specified for the `mode` parameter). In the
case where the user wants to flag multiple records of a grouping, for
example records that all happen on the same visit and time, the argument
`flag_all` can be set to `TRUE`. Otherwise, `new_var` is set to `NA`.
Thus, the direction of "worst" is considered fixed for all parameters in
the dataset depending on the `order` and the `mode`, i.e. for every
parameter the first or last record will be flagged across the whole
dataset.

## See also

General Derivation Functions for all ADaMs that returns variable
appended to dataset:
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_joined_exist_flag.md),
[`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_ef_msrc.md),
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_merged_exist_flag.md),
[`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_obs_number.md),
[`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_relative_flag.md),
[`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_cat.md),
[`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_computed.md),
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined.md),
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_joined_summary.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged.md),
[`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_lookup.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_merged_summary.md),
[`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_transposed.md)

## Examples

### Data setup

The following examples use the `ADVS` and `ADAE` datasets below as a
basis.

    library(tibble, warn.conflicts = FALSE)
    library(lubridate, warn.conflicts = FALSE)
    library(dplyr, warn.conflicts = FALSE)

    advs <- tribble(
      ~USUBJID, ~PARAMCD,    ~AVISIT,          ~ADT, ~AVAL,
      "1015",   "TEMP",   "BASELINE",  "2021-04-27",  38.0,
      "1015",   "TEMP",   "BASELINE",  "2021-04-25",  39.0,
      "1015",   "TEMP",   "WEEK 2",    "2021-05-10",  37.5,
      "1015",   "WEIGHT", "SCREENING", "2021-04-19",  81.2,
      "1015",   "WEIGHT", "BASELINE",  "2021-04-25",  82.7,
      "1015",   "WEIGHT", "BASELINE",  "2021-04-27",  84.0,
      "1015",   "WEIGHT", "WEEK 2",    "2021-05-09",  82.5,
      "1023",   "TEMP",   "SCREENING", "2021-04-27",  38.0,
      "1023",   "TEMP",   "BASELINE",  "2021-04-28",  37.5,
      "1023",   "TEMP",   "BASELINE",  "2021-04-29",  37.5,
      "1023",   "TEMP",   "WEEK 1",    "2021-05-03",  37.0,
      "1023",   "WEIGHT", "SCREENING", "2021-04-27",  69.6,
      "1023",   "WEIGHT", "BASELINE",  "2021-04-29",  67.2,
      "1023",   "WEIGHT", "WEEK 1",    "2021-05-02",  65.9
    ) %>%
    mutate(
      STUDYID = "AB123",
      ADT = ymd(ADT)
    )

    adae <- tribble(
      ~USUBJID,         ~AEBODSYS,    ~AEDECOD,   ~AESEV, ~AESTDY, ~AESEQ,
      "1015", "GENERAL DISORDERS",  "ERYTHEMA",   "MILD",       2,      1,
      "1015", "GENERAL DISORDERS",  "PRURITUS",   "MILD",       2,      2,
      "1015",      "GI DISORDERS", "DIARRHOEA",   "MILD",       8,      3,
      "1023", "CARDIAC DISORDERS",  "AV BLOCK",   "MILD",      22,      4,
      "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       3,      1,
      "1023",    "SKIN DISORDERS",  "ERYTHEMA", "SEVERE",       5,      2,
      "1023",    "SKIN DISORDERS",  "ERYTHEMA",   "MILD",       8,      3
    ) %>%
    mutate(STUDYID = "AB123")

### Flagging the first/last observation within a by group (`order`, `mode`)

A new variable is added for each subject to flag the last observation
within a by group. Within each by group (specified by `by_vars`), the
`order = exprs(ADT)` argument specifies we wish to sort the records by
analysis date and then select the last one (`mode = "last"`). The name
of the new variable is passed through the `new_var = LASTFL` call.

    advs %>%
      derive_var_extreme_flag(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(ADT),
        new_var = LASTFL,
        mode = "last",
      ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL LASTFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   <NA>
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   <NA>
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 Y
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 <NA>
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 <NA>
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   <NA>
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 Y
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   <NA>
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 <NA>
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 <NA>
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   Y
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 <NA>
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 <NA>
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 Y     

Note here that a similar `FIRSTFL` variable could instead be derived
simply by switching to `mode = "first"`. Alternatively, we could make
use of
[`desc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-desc.md)
within the sorting specified by `order`:

    advs %>%
      derive_var_extreme_flag(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(desc(ADT)),
        new_var = FIRSTFL,
        mode = "last",
      ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL FIRSTFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   Y
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   <NA>
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 <NA>
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 Y
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 <NA>
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   <NA>
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 <NA>
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   Y
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 <NA>
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 <NA>
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   <NA>
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 Y
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 <NA>
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 <NA>   

### Modifying the flag values (`true_value`, `false_value`)

The previous example is now enhanced with custom values for the flag
entries. Records which are flagged are filled with the contents of
`true_value` and those which are not are filled with the contents of
`false_value`. Note that these are normally preset to `"Y"` and `NA`,
which is why they were not specified in the example above.

    advs %>%
      derive_var_extreme_flag(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(ADT),
        new_var = LASTFL,
        mode = "last",
        true_value = "Yes",
        false_value = "No",
      ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL LASTFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   No
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   No
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 Yes
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 No
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 No
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   No
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 Yes
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   No
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 No
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 No
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   Yes
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 No
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 No
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 Yes   

### Creating temporary variables for sorting (`check_type`)

In this example we wish to flag the first occurrence of the most severe
AE within each subject. To ensure correct sorting of the severity
values, `AESEV` must be pre-processed into a numeric variable
`TEMP_AESEVN` which can then be passed inside `order`. Once again, to
ensure we only flag the *first* occurrence, we specify `AESTDY` and
`AESEQ` inside `order` as well.

    adae %>%
      mutate(
        TEMP_AESEVN =
          as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
      ) %>%
      derive_var_extreme_flag(
        new_var = AOCCIFL,
        by_vars = exprs(STUDYID, USUBJID),
        order = exprs(TEMP_AESEVN, AESTDY, AESEQ),
        mode = "first",
        check_type = "warning"
      ) %>%
      arrange(STUDYID, USUBJID, AESTDY, AESEQ) %>%
      select(STUDYID, USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
    #> # A tibble: 7 × 7
    #>   STUDYID USUBJID AEDECOD   AESEV  AESTDY AESEQ AOCCIFL
    #>   <chr>   <chr>   <chr>     <chr>   <dbl> <dbl> <chr>
    #> 1 AB123   1015    ERYTHEMA  MILD        2     1 Y
    #> 2 AB123   1015    PRURITUS  MILD        2     2 <NA>
    #> 3 AB123   1015    DIARRHOEA MILD        8     3 <NA>
    #> 4 AB123   1023    ERYTHEMA  MILD        3     1 <NA>
    #> 5 AB123   1023    ERYTHEMA  SEVERE      5     2 Y
    #> 6 AB123   1023    ERYTHEMA  MILD        8     3 <NA>
    #> 7 AB123   1023    AV BLOCK  MILD       22     4 <NA>   

Note here that the presence of `AESEQ` as a sorting variable inside the
`order` argument ensures that the combination of `by_vars` and `order`
indexes unique records in the dataset. If this had been omitted, the
choice of `check_type = "warning"` would have ensured that
`derive_var_extreme_flag()` would throw a warning due to perceived
duplicate records (in this case, the first two AEs for subject
`"1015"`). If no sorting variables exist, or if these duplicates are
acceptable, then the user can silence the warning with
`check_type = "none"`. Alternatively, the warning can be upgraded to an
error with `check_type = "error"`.

### Flagging all records if multiple are identified (`flag_all`)

Revisiting the above example, if we instead wish to flag *all* AEs of
the highest severity occurring on the earliest date, then we can use
`flag_all = TRUE`. Note that we now also omit `AESEQ` from the `order`
argument because we do not need to differentiate between two AEs
occurring on the same day (e.g. for subject `"1015"`) as they are both
flagged.

    adae %>%
      mutate(
        TEMP_AESEVN =
          as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))
      ) %>%
      derive_var_extreme_flag(
        new_var = AOCCIFL,
        by_vars = exprs(STUDYID, USUBJID),
        order = exprs(TEMP_AESEVN, AESTDY),
        mode = "first",
        flag_all = TRUE
      ) %>%
      arrange(STUDYID, USUBJID, AESTDY, AESEQ) %>%
      select(STUDYID, USUBJID, AEDECOD, AESEV, AESTDY, AESEQ, AOCCIFL)
    #> # A tibble: 7 × 7
    #>   STUDYID USUBJID AEDECOD   AESEV  AESTDY AESEQ AOCCIFL
    #>   <chr>   <chr>   <chr>     <chr>   <dbl> <dbl> <chr>
    #> 1 AB123   1015    ERYTHEMA  MILD        2     1 Y
    #> 2 AB123   1015    PRURITUS  MILD        2     2 Y
    #> 3 AB123   1015    DIARRHOEA MILD        8     3 <NA>
    #> 4 AB123   1023    ERYTHEMA  MILD        3     1 <NA>
    #> 5 AB123   1023    ERYTHEMA  SEVERE      5     2 Y
    #> 6 AB123   1023    ERYTHEMA  MILD        8     3 <NA>
    #> 7 AB123   1023    AV BLOCK  MILD       22     4 <NA>   

### Deriving a baseline flag

`derive_var_extreme_flag()` is very often used to derive the baseline
flag `ABLFL`, so the following section contains various examples of this
in action for the `ADVS` dataset. Note that for these derivations it is
often convenient to leverage higher order functions such as
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
and
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/slice_derivation.md).
Please read the [Higher Order
Functions](https://pharmaverse.github.io/admiral/articles/higher_order.html)
vignette, as well as their specific reference pages, to learn more.

To set the baseline flag for the last observation among those where
`AVISIT = "BASELINE"`, we can use a similar call to the examples above
but wrapping inside of
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
and making use of the `filter` argument.

    restrict_derivation(
      advs,
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(USUBJID, PARAMCD),
        order = exprs(ADT),
        new_var = ABLFL,
        mode = "last"
      ),
      filter = AVISIT == "BASELINE"
    ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL ABLFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   <NA>
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   Y
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 <NA>
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 <NA>
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 <NA>
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   Y
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 <NA>
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   <NA>
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 <NA>
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 Y
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   <NA>
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 <NA>
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 Y
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 <NA> 

Alternatively, to set baseline as the lowest observation among those
where `AVISIT = "BASELINE"` (selecting the latest if there are multiple)
we can modify the `order` argument, ensuring to sort by descending
`AVAL` before `ADT`. Note here the synergy between
[`desc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-desc.md)
and `mode`, because `mode = "last"` applies to both the ordering
variables `AVAL` and `ADT` and so we need to reverse only the ordering
of the former to ensure that the lowest value is selected but also that
the latest one among multiple is preferred. This is relevant for subject
`"1023"`'s temperature records.

    restrict_derivation(
      advs,
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(USUBJID, PARAMCD),
        order = exprs(desc(AVAL), ADT),
        new_var = ABLFL,
        mode = "last"
      ),
      filter = AVISIT == "BASELINE"
    ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL ABLFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   <NA>
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   Y
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 <NA>
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 <NA>
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 Y
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   <NA>
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 <NA>
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   <NA>
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 <NA>
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 Y
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   <NA>
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 <NA>
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 Y
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 <NA> 

In practice, baseline-setting may vary on a parameter by parameter
basis, in which case
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/slice_derivation.md)
could be used in place of
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md).
In the example below, we set the baseline flag as follows: for
temperature records, as the lowest value recorded at a baseline visit;
for weight records, as the highest value recorded at a baseline visit.
In both cases, we again select the latest observation if there are
multiple.

    slice_derivation(
      advs,
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(USUBJID, PARAMCD),
        mode = "last",
        new_var = ABLFL,
      ),
      derivation_slice(
        filter = AVISIT == "BASELINE" & PARAMCD == "TEMP",
        args = params(order = exprs(desc(AVAL), ADT))
      ),
      derivation_slice(
        filter = AVISIT == "BASELINE" & PARAMCD == "WEIGHT",
        args = params(order = exprs(AVAL, ADT))
      )
    ) %>%
      arrange(STUDYID, USUBJID, PARAMCD, ADT) %>%
      select(STUDYID, everything())
    #> # A tibble: 14 × 7
    #>    STUDYID USUBJID PARAMCD AVISIT    ADT         AVAL ABLFL
    #>    <chr>   <chr>   <chr>   <chr>     <date>     <dbl> <chr>
    #>  1 AB123   1015    TEMP    BASELINE  2021-04-25  39   <NA>
    #>  2 AB123   1015    TEMP    BASELINE  2021-04-27  38   Y
    #>  3 AB123   1015    TEMP    WEEK 2    2021-05-10  37.5 <NA>
    #>  4 AB123   1015    WEIGHT  SCREENING 2021-04-19  81.2 <NA>
    #>  5 AB123   1015    WEIGHT  BASELINE  2021-04-25  82.7 <NA>
    #>  6 AB123   1015    WEIGHT  BASELINE  2021-04-27  84   Y
    #>  7 AB123   1015    WEIGHT  WEEK 2    2021-05-09  82.5 <NA>
    #>  8 AB123   1023    TEMP    SCREENING 2021-04-27  38   <NA>
    #>  9 AB123   1023    TEMP    BASELINE  2021-04-28  37.5 <NA>
    #> 10 AB123   1023    TEMP    BASELINE  2021-04-29  37.5 Y
    #> 11 AB123   1023    TEMP    WEEK 1    2021-05-03  37   <NA>
    #> 12 AB123   1023    WEIGHT  SCREENING 2021-04-27  69.6 <NA>
    #> 13 AB123   1023    WEIGHT  BASELINE  2021-04-29  67.2 Y
    #> 14 AB123   1023    WEIGHT  WEEK 1    2021-05-02  65.9 <NA> 
