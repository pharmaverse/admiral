# Generic Derivations

## Introduction

This vignette provides an overview of the most important admiral
functions, the generic derivations. They do not derive a specific
variable or parameter in a specific dataset but can be used for many
different derivations in many different datasets. The most important
concepts and some examples are presented here. For full details and more
examples see the documentation of the functions.

### Required Packages

The examples in this vignette require the following packages.

``` r
library(admiral)
library(tibble)
```

## Characterization of Derivations

The generic functions can be characterized by the following three
properties:

- What should be added? There are two options:
  - *Variables*
  - *Records/Parameters*
- What is the source data? There are two options:
  - A *single* source dataset
  - *Multiple* source datasets
- Which method should be used? There are three options:
  - *Selection*: The new values are derived by selecting records, e.g.,
    the baseline records, the last exposure record, …
  - *Summary*: The new values are derived by summarizing values, e.g.,
    sum, average, concatenation, …
  - *Computation*: The new values are derived by a computation with more
    than one value as input, e.g., deriving BMI or BSA from height and
    weight

### Overview of Derivations

Using the three properties makes it easy to find the appropriate
function for a particular derivation. The following interactive table
lists all generic functions and their properties.

## Source Data

Most derivation functions expect a single source dataset. For some
multiple source datasets can be specified. In both cases the way how to
specify the source datasets is the same across all generic functions.

### Single Source Dataset

For functions expecting a single source dataset the data is provided by
the `dataset_add` argument. This is a mandatory argument. The data
provided by the `dataset` argument is not used[¹](#fn1).

If the `dataset_add` argument is not provided, the data from `dataset`
is used
([`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_flag.md)).

### Multiple Source Datasets

For functions expecting multiple source datasets the data is provided by
the `source_datasets` argument. The datasets are referred to by the
`dataset_name` element of the source objects.

For example, consider the derivation of a response parameter. The three
possible responses are defined by
[`event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event.md)
objects. These objects define the events but do not include any data.
Instead the `dataset_name` field is set to a (character) id. This id is
used in the `source_datasets` argument of the derivation function to
link the data to the events. I.e., for the first two events
(`complete_response` and `partial_response`) the dataset `adrs_ovr` is
used while for the last event the dataset `myadsl` is used.

    complete_response <- event(
      description = "Define complete response",
      dataset_name = "ovr",
      condition = AVALC == "CR",
      set_values_to = exprs(AVALC = "COMPLETE RESPONSE")
    )

    partial_response <- event(
      description = "Define partial response",
      dataset_name = "ovr",
      condition = AVALC == "PR",
      set_values_to = exprs(AVALC = "PARTIAL RESPONSE")
    )

    no_response <- event(
      description = "Define no response for all patients in adsl",
      dataset_name = "adsl",
      condition = TRUE,
      set_values_to = exprs(AVALC = "NO RESPONSE")
    )

    derive_extreme_event(
      ...
      events = list(complete_response, partial_response, no_response),
      source_datasets = list(ovr = adrs_ovr, adsl = myadsl),
      ...
    )

This allows to define the source objects independent of the data. I.e.,
the same source object can be used for different source datasets. For
example, the parameter above could be derived for data from a second
reporter by just changing `source_dataset`:

    derive_extreme_event(
      ...
      events = list(complete_response, partial_response, no_response),
      source_datasets = list(ovr = adrs_ovr_reporter2, adsl = myadsl),
      ...
    )

For some source objects the `dataset_name` element is optional, e.g.,
[`event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/event.md).
If it is not specified, the input dataset (`dataset`) is used.

## Methods

The generic derivations use three different methods for deriving the
values of the new variables or records. Usually the derivation applies
the method several times, once for each group in the input or source
data. The groups are defined by the `by_vars` argument, e.g.,
`by_vars = exprs(USUBJID)` for “each subject” or
`by_vars = exprs(USUBJID, PARAMCD)` for “each subject and parameter”.

### Selection

The most common method is the selection method. It selects a record from
the source dataset(s) and adds information from the selected record to
the input dataset. This information could be just a flag indicating that
a record exists, one or more variables from the selected records, or new
variables created from the selected record.

#### Options for Selection

In the simplest case the record is selected by a condition. The
condition is specified by the `filter_add` or `filter` argument. In the
following example baseline weight is added to `ADSL`.

``` r
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
)

advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~ABLFL, ~AVAL, ~AVALU,
  "1",      "WEIGHT", "BASELINE", "Y",     58.7, "kg",
  "1",      "HEIGHT", "BASELINE", "Y",    169.2, "cm",
  "1",      "WEIGHT", "WEEK 3",   NA,      59.3, "kg",
  "2",      "WEIGHT", "BASELINE", "Y",     72.5, "kg",
  "2",      "WEIGHT", "WEKK 3",   NA,      71.9, "kg",
)

derive_vars_merged(
  adsl,
  dataset_add = advs,
  by_vars = exprs(USUBJID),
  filter_add = PARAMCD == "WEIGHT" & ABLFL == "Y",
  new_vars = exprs(WGTBL = AVAL)
)
#> # A tibble: 3 × 2
#>   USUBJID WGTBL
#>   <chr>   <dbl>
#> 1 1        58.7
#> 2 2        72.5
#> 3 3        NA
```

Sometimes it is not possible to select the record of interest by a
condition, e.g., if the first, last, best, worst, lowest, highest, …
value should be derived. In this case the `mode` and `order` argument
can be specified to select the first or last record with respect to the
variables specified for `order`. Below the day of the last valid dose is
added to `ADSL`.

``` r
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
)

ex <- tribble(
  ~USUBJID, ~EXSTDY, ~EXDOSE,
  "1",            1,      50,
  "1",            7,      70,
  "1",           14,       0,
  "2",            1,      75,
  "2",            9,      70
)

derive_vars_merged(
  adsl,
  dataset_add = ex,
  by_vars = exprs(USUBJID),
  filter_add = EXDOSE > 0,
  order = exprs(EXSTDY),
  mode = "last",
  new_vars = exprs(TRTEDY = EXSTDY)
)
#> # A tibble: 3 × 2
#>   USUBJID TRTEDY
#>   <chr>    <dbl>
#> 1 1            7
#> 2 2            9
#> 3 3           NA
```

It is also possible to select the record based on records of the input
*and* the source dataset. For this type of selection
[admiral](https://pharmaverse.github.io/admiral/) provides the functions
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined.md),
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md),
and
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md).
They provide the `filter_join` argument which accepts conditions with
variables from both the input dataset (`dataset`) and the additional
dataset (`dataset_add`). As an example consider deriving the day and
dose of the last study treatment before an adverse event:

``` r
adae <- tribble(
  ~USUBJID, ~ASTDY, ~AESEQ,
  "1",           3,      1,
  "1",           3,      2,
  "1",          15,      3
)

ex <- tribble(
  ~USUBJID, ~EXSTDY, ~EXDOSE,
  "1",            1,      50,
  "1",            7,      70,
  "1",           14,       0,
  "2",            1,      75,
  "2",            9,      70
)

derive_vars_joined(
  adae,
  dataset_add = ex,
  by_vars = exprs(USUBJID),
  filter_add = EXDOSE > 0,
  filter_join = EXSTDY <= ASTDY,
  join_type = "all",
  order = exprs(EXSTDY),
  mode = "last",
  new_vars = exprs(LSTDOSDY = EXSTDY, LASTDOS = EXDOSE)
)
#> # A tibble: 3 × 5
#>   USUBJID ASTDY AESEQ LSTDOSDY LASTDOS
#>   <chr>   <dbl> <dbl>    <dbl>   <dbl>
#> 1 1           3     1        1      50
#> 2 1           3     2        1      50
#> 3 1          15     3        7      70
```

The `filter_join` condition is applied on a temporary dataset created by
left joining the input dataset and the additional dataset (restricted by
`filter_add`):

    #> # A tibble: 6 × 5
    #>   USUBJID ASTDY AESEQ EXSTDY EXDOSE
    #>   <chr>   <dbl> <dbl>  <dbl>  <dbl>
    #> 1 1           3     1      1     50
    #> 2 1           3     1      7     70
    #> 3 1           3     2      1     50
    #> 4 1           3     2      7     70
    #> 5 1          15     3      1     50
    #> 6 1          15     3      7     70

The “joined” function can also be used when the condition for selecting
depends on previous or subsequent records in the dataset. In this case
the same dataset is specified for `dataset` and `dataest_add`. Consider
the following example where `"HIGH"` results should be flagged which are
confirmed by a second `"HIGH"` result at least ten days later.

``` r
adlb <- tribble(
  ~USUBJID, ~PARAMCD, ~ADY, ~ANRIND,
  "1",      "AST",       1, "HIGH",
  "1",      "AST",       7, "HIGH",
  "1",      "AST",      14, "NORMAL",
  "1",      "ALT",       1, "HIGH",
  "1",      "ALT",       7, "NORMAL",
  "1",      "ALT",      14, "HIGH",
  "2",      "AST",       1, "HIGH",
  "2",      "AST",      15, "HIGH",
  "2",      "AST",      22, "NORMAL",
  "2",      "ALT",       1, "HIGH"
)

derive_var_joined_exist_flag(
  adlb,
  dataset_add = adlb,
  by_vars = exprs(USUBJID, PARAMCD),
  order = exprs(ADY),
  join_vars = exprs(ADY, ANRIND),
  join_type = "after",
  filter_join = ANRIND == "HIGH" & ANRIND.join == "HIGH" & ADY.join > ADY + 10,
  new_var = HICONFFL
)
#> # A tibble: 10 × 5
#>    USUBJID PARAMCD   ADY ANRIND HICONFFL
#>    <chr>   <chr>   <dbl> <chr>  <chr>   
#>  1 1       AST         1 HIGH   NA      
#>  2 1       AST         7 HIGH   NA      
#>  3 1       AST        14 NORMAL NA      
#>  4 1       ALT         1 HIGH   Y       
#>  5 1       ALT         7 NORMAL NA      
#>  6 1       ALT        14 HIGH   NA      
#>  7 2       AST         1 HIGH   Y       
#>  8 2       AST        15 HIGH   NA      
#>  9 2       AST        22 NORMAL NA      
#> 10 2       ALT         1 HIGH   NA
```

The `join_type` argument is set to `"after"` to restrict the joined
records to subsequent results.

As the same variables are included in `dataset` and `dataset_add`, those
from `dataset_add` are renamed by adding the suffix “.join”. The
variables from `dataset_add` which are used in `filter_join` must be
specified for `join_vars`. So the temporary dataset for applying
`filter_join` is:

    #> # A tibble: 9 × 6
    #>   USUBJID PARAMCD   ADY ANRIND ADY.join ANRIND.join
    #>   <chr>   <chr>   <dbl> <chr>     <dbl> <chr>      
    #> 1 1       ALT         1 HIGH          7 NORMAL     
    #> 2 1       ALT         1 HIGH         14 HIGH       
    #> 3 1       ALT         7 NORMAL       14 HIGH       
    #> 4 1       AST         1 HIGH          7 HIGH       
    #> 5 1       AST         1 HIGH         14 NORMAL     
    #> 6 1       AST         7 HIGH         14 NORMAL     
    #> 7 2       AST         1 HIGH         15 HIGH       
    #> 8 2       AST         1 HIGH         22 NORMAL     
    #> 9 2       AST        15 HIGH         22 NORMAL

It is possible to use summary functions like
[`all()`](https://rdrr.io/r/base/all.html) or
[`any()`](https://rdrr.io/r/base/any.html) in `filter_join`. Assume that
in the previous example records should be flagged only if all results
between the flagged record and the confirmation record were `"HIGH"`.
This can be achieved by specifying the `first_cond_upper` argument and
set it to the condition for confirmation.

``` r
derive_var_joined_exist_flag(
  adlb,
  dataset_add = adlb,
  by_vars = exprs(USUBJID, PARAMCD),
  order = exprs(ADY),
  join_vars = exprs(ADY, ANRIND),
  join_type = "after",
  first_cond_upper = ANRIND.join == "HIGH" & ADY.join > ADY + 10,
  filter_join = ANRIND == "HIGH" & all(ANRIND.join == "HIGH"),
  new_var = HICONFFL
)
#> # A tibble: 10 × 5
#>    USUBJID PARAMCD   ADY ANRIND HICONFFL
#>    <chr>   <chr>   <dbl> <chr>  <chr>   
#>  1 1       AST         1 HIGH   NA      
#>  2 1       AST         7 HIGH   NA      
#>  3 1       AST        14 NORMAL NA      
#>  4 1       ALT         1 HIGH   NA      
#>  5 1       ALT         7 NORMAL NA      
#>  6 1       ALT        14 HIGH   NA      
#>  7 2       AST         1 HIGH   Y       
#>  8 2       AST        15 HIGH   NA      
#>  9 2       AST        22 NORMAL NA      
#> 10 2       ALT         1 HIGH   NA
```

If the `first_cond_upper` argument is specified the records in the
joined dataset are restricted up to the first records where the
condition is fulfilled:

    #> # A tibble: 3 × 6
    #>   USUBJID PARAMCD   ADY ANRIND ADY.join ANRIND.join
    #>   <chr>   <chr>   <dbl> <chr>     <dbl> <chr>      
    #> 1 1       ALT         1 HIGH          7 NORMAL     
    #> 2 1       ALT         1 HIGH         14 HIGH       
    #> 3 2       AST         1 HIGH         15 HIGH

Thereafter `filter_join` is applied to the restricted joined dataset.
I.e., the [`all()`](https://rdrr.io/r/base/all.html) function considers
only the results up to the confirmation records and ignores subsequent
results.

**Note:** In principle, we actually could achieve every result from
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md)
and
[`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_exist_flag.md)
by using
[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined.md)
or
[`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md)
respectively. However, the “joined” functions require much more
resources (time and memory), hence it is recommended to use them only if
it is really required, i.e., the condition for selecting records depends
on variables from *both* datasets.

#### Sort Order

The [admiral](https://pharmaverse.github.io/admiral/) functions use
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)
for sorting, i.e., `NA`s are always sorted to the end (regardless
whether
[`desc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-desc.md)
is used or not).

Consider for example the following derivation of a last visit flag. The
record with `AVISITN == NA` is flagged because `NA` is sorted to the
end.

``` r
advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVISITN, ~AVAL,
  "1",      "WEIGHT",       NA,  62.1,
  "1",      "WEIGHT",        1,  62.3,
  "1",      "WEIGHT",        2,  62.5,
  "1",      "WEIGHT",        3,  62.4
)

derive_var_extreme_flag(
  advs,
  by_vars = exprs(USUBJID, PARAMCD),
  order = exprs(AVISITN),
  mode = "last",
  new_var = LSTVISFL
)
#> # A tibble: 4 × 5
#>   USUBJID PARAMCD AVISITN  AVAL LSTVISFL
#>   <chr>   <chr>     <dbl> <dbl> <chr>   
#> 1 1       WEIGHT        1  62.3 NA      
#> 2 1       WEIGHT        2  62.5 NA      
#> 3 1       WEIGHT        3  62.4 NA      
#> 4 1       WEIGHT       NA  62.1 Y
```

The `order` argument accepts expressions. This allows to specify how
`NA`s should be handled. For example, the following sorts the `NA` to
the start. Thus the `AVISITN == 3` record is flagged.

``` r
derive_var_extreme_flag(
  advs,
  by_vars = exprs(USUBJID, PARAMCD),
  order = exprs(if_else(is.na(AVISITN), -Inf, AVISITN)),
  mode = "last",
  new_var = LSTVISFL
)
#> # A tibble: 4 × 5
#>   USUBJID PARAMCD AVISITN  AVAL LSTVISFL
#>   <chr>   <chr>     <dbl> <dbl> <chr>   
#> 1 1       WEIGHT       NA  62.1 NA      
#> 2 1       WEIGHT        1  62.3 NA      
#> 3 1       WEIGHT        2  62.5 NA      
#> 4 1       WEIGHT        3  62.4 Y
```

The same can achieved with the following, which also works for character
variables.

``` r
derive_var_extreme_flag(
  advs,
  by_vars = exprs(USUBJID, PARAMCD),
  order = exprs(!is.na(AVISITN), AVISITN),
  mode = "last",
  new_var = LSTVISFL
)
#> # A tibble: 4 × 5
#>   USUBJID PARAMCD AVISITN  AVAL LSTVISFL
#>   <chr>   <chr>     <dbl> <dbl> <chr>   
#> 1 1       WEIGHT       NA  62.1 NA      
#> 2 1       WEIGHT        1  62.3 NA      
#> 3 1       WEIGHT        2  62.5 NA      
#> 4 1       WEIGHT        3  62.4 Y
```

#### New Values

How the (new) variables are set depends on whether variables, a flag, or
records are added by the derivation.

- If only a flag needs to be added, the flag functions
  ([`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_exist_flag.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_joined_exist_flag.md))
  can be used. The name of the new variable is specified by the
  `new_var` argument and the values of the flag by `true_value` and
  `false_value`.
- If new variables from the selected record needs to be added, the name
  of the new variables and their values are specified by the `new_vars`
  argument. If in addition a flag should be added, the `exist_flag`
  argument and the `true_value` and `false_value` argument can be used.
- If new records are added, the variables and their values are defined
  by the `set_values_to` argument.

### Summary

If the new values should be derived by summarizing values, e.g., sum,
average, concatenation, …, the functions
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md),
[`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged_summary.md),
or
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined_summary.md)
can be used. For example, adding an average dose parameter to `ADEX` can
be done by the following:

``` r
adex <- tribble(
  ~USUBJID, ~ASTDY, ~AVAL, ~PARAMCD,
  "1",           1,    50, "DOSE",
  "1",           7,    70, "DOSE",
  "1",          14,     0, "DOSE",
  "2",           1,    75, "DOSE",
  "2",           9,    70, "DOSE"
)

derive_summary_records(
  adex,
  dataset_add = adex,
  filter_add = AVAL > 0,
  by_vars = exprs(USUBJID),
  set_values_to = exprs(
    AVAL = mean(AVAL),
    PARAMCD = "AVERAGE DOSE"
  )
)
#> # A tibble: 7 × 4
#>   USUBJID ASTDY  AVAL PARAMCD     
#>   <chr>   <dbl> <dbl> <chr>       
#> 1 1           1  50   DOSE        
#> 2 1           7  70   DOSE        
#> 3 1          14   0   DOSE        
#> 4 2           1  75   DOSE        
#> 5 2           9  70   DOSE        
#> 6 1          NA  60   AVERAGE DOSE
#> 7 2          NA  72.5 AVERAGE DOSE
```

The summary function, the source variable, and the new variable can be
specified by the `set_values_to` argument. Variables which are not
specified by `set_values_to` or `by_vars` are set to `NA` for the new
records.

If the average dose should be added as a variable to `ADSL`, consider
the following:

``` r
adsl <- tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
)

derive_vars_merged_summary(
  adsl,
  dataset_add = adex,
  filter_add = AVAL > 0,
  by_vars = exprs(USUBJID),
  new_vars = exprs(
    AVERDOSE = mean(AVAL)
  ),
  missing_values = exprs(AVERDOSE = 0)
)
#> # A tibble: 3 × 2
#>   USUBJID AVERDOSE
#>   <chr>      <dbl>
#> 1 1           60  
#> 2 2           72.5
#> 3 3            0
```

Here the summary function, the source variable, and the new variable are
specified by the `new_vars` argument. For subjects without exposure
observations the value of the new variable can be defined by the
`missing_values` argument.

If the selection of records to summarize depends on the records of both
the input dataset and the additional dataset, the
[`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_joined_summary.md)
function can be used. For example, the following adds the cumulative
dose up to the event as a variable (`CUMDOSA`) to the input dataset:

``` r
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
  new_vars = exprs(CUMDOSA = sum(AVAL, na.rm = TRUE)),
  missing_value = exprs(CUMDOSA = 0)
)
#> # A tibble: 6 × 4
#>   USUBJID   ADY AEDECOD    CUMDOSA
#>   <chr>   <dbl> <chr>        <dbl>
#> 1 1           2 Fatigue         10
#> 2 1           9 Influenza       30
#> 3 1          15 Theft           40
#> 4 1          15 Fatigue         40
#> 5 2           4 Parasomnia       0
#> 6 3           2 Truancy          0
```

### Computed

If the new values should be computed from different parameters of the
source dataset,
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md)
can be used. The computed value can be specified by the `set_values_to`
argument. The values of the source variable for a parameter can be
referred to by temporary variables of the form
`<variable name>.<parameter name>`.

``` r
advs <- tribble(
  ~USUBJID, ~AVISIT,    ~PARAMCD, ~AVAL, ~AVALU,
  "1",      "BASELINE", "WEIGHT",  32.6, "kg",
  "1",      "BASELINE", "HEIGHT", 155.4, "cm",
  "1",      "MONTH 6",  "WEIGHT",  33.2, "kg",
  "1",      "MONTH 6",  "HEIGHT", 155.8, "cm",
  "2",      "BASELINE", "WEIGHT",  44.2, "kg",
  "2",      "BASELINE", "HEIGHT", 145.3, "cm",
  "2",      "MONTH 6",  "WEIGHT",  42.0, "kg",
  "2",      "MONTH 6",  "HEIGHT", 146.4, "cm"
)

derive_param_computed(
  advs,
  by_vars = exprs(USUBJID, AVISIT),
  parameters = c("WEIGHT", "HEIGHT"),
  set_values_to = exprs(
    AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
    PARAMCD = "BMI",
    AVALU = "kg/m^2"
  )
)
#> # A tibble: 12 × 5
#>    USUBJID AVISIT   PARAMCD  AVAL AVALU 
#>    <chr>   <chr>    <chr>   <dbl> <chr> 
#>  1 1       BASELINE WEIGHT   32.6 kg    
#>  2 1       BASELINE HEIGHT  155.  cm    
#>  3 1       MONTH 6  WEIGHT   33.2 kg    
#>  4 1       MONTH 6  HEIGHT  156.  cm    
#>  5 2       BASELINE WEIGHT   44.2 kg    
#>  6 2       BASELINE HEIGHT  145.  cm    
#>  7 2       MONTH 6  WEIGHT   42   kg    
#>  8 2       MONTH 6  HEIGHT  146.  cm    
#>  9 1       BASELINE BMI      13.5 kg/m^2
#> 10 1       MONTH 6  BMI      13.7 kg/m^2
#> 11 2       BASELINE BMI      20.9 kg/m^2
#> 12 2       MONTH 6  BMI      19.6 kg/m^2
```

For common computations like BMI
[admiral](https://pharmaverse.github.io/admiral/) offers [computation
functions](https://pharmaverse.github.io/admiral/cran-release/reference/#computation-functions-for-vectors).
In the previous example
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_bmi.md)
could be used instead of the formula for BMI:

``` r
derive_param_computed(
  advs,
  by_vars = exprs(USUBJID, AVISIT),
  parameters = c("WEIGHT", "HEIGHT"),
  set_values_to = exprs(
    AVAL = compute_bmi(weight = AVAL.WEIGHT, height = AVAL.HEIGHT),
    PARAMCD = "BMI",
    AVALU = "kg/m^2"
  )
)
#> # A tibble: 12 × 5
#>    USUBJID AVISIT   PARAMCD  AVAL AVALU 
#>    <chr>   <chr>    <chr>   <dbl> <chr> 
#>  1 1       BASELINE WEIGHT   32.6 kg    
#>  2 1       BASELINE HEIGHT  155.  cm    
#>  3 1       MONTH 6  WEIGHT   33.2 kg    
#>  4 1       MONTH 6  HEIGHT  156.  cm    
#>  5 2       BASELINE WEIGHT   44.2 kg    
#>  6 2       BASELINE HEIGHT  145.  cm    
#>  7 2       MONTH 6  WEIGHT   42   kg    
#>  8 2       MONTH 6  HEIGHT  146.  cm    
#>  9 1       BASELINE BMI      13.5 kg/m^2
#> 10 1       MONTH 6  BMI      13.7 kg/m^2
#> 11 2       BASELINE BMI      20.9 kg/m^2
#> 12 2       MONTH 6  BMI      19.6 kg/m^2
```

------------------------------------------------------------------------

1.  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_computed.md)
    is an exception. It uses the data from both arguments.
