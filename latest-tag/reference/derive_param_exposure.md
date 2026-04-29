# Add an Aggregated Parameter and Derive the Associated Start and End Dates

Add a record computed from the aggregated analysis value of another
parameter and compute the start (`ASTDT(M)`)and end date (`AENDT(M)`) as
the minimum and maximum date by `by_vars`.

## Usage

``` r
derive_param_exposure(
  dataset = NULL,
  dataset_add,
  by_vars,
  input_code,
  filter_add = NULL,
  set_values_to = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset.

  Default value

  :   `NULL`

- dataset_add:

  Additional dataset

  The variables specified for `by_vars`, `analysis_var`, `PARAMCD`,
  alongside either `ASTDTM` and `AENDTM` or `ASTDT` and `AENDT` are also
  expected. Observations from the specified dataset are going to be used
  to calculate and added as new records to the input dataset
  (`dataset`).

  Default value

  :   none

- by_vars:

  Grouping variables

  For each group defined by `by_vars` an observation is added to the
  output dataset. Only variables specified in `by_vars` will be
  populated in the newly created records.

  Default value

  :   none

- input_code:

  Required parameter code

  The observations where `PARAMCD` equals the specified value are
  considered to compute the summary record.

  Permitted values

  :   A character of `PARAMCD` value

  Default value

  :   none

- filter_add:

  Filter condition as logical expression to apply during summary
  calculation. By default, filtering expressions are computed within
  `by_vars` as this will help when an aggregating, lagging, or ranking
  function is involved.

  For example,

  - `filter_add = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all
    `AVAL` values greater than mean of `AVAL` with in `by_vars`.

  - `filter_add = (dplyr::n() > 2)` will filter n count of `by_vars`
    greater than 2.

  Default value

  :   `NULL`

- set_values_to:

  Variable-value pairs

  Set a list of variables to some specified value for the new
  observation(s)

  - LHS refer to a variable. It is expected that at least `PARAMCD` is
    defined.

  - RHS refers to the values to set to the variable. This can be a
    string, a symbol, a numeric value, `NA`, or an expression. (e.g.
    `exprs(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")`).

  Permitted values

  :   List of variable-value pairs

  Default value

  :   `NULL`

## Value

The input dataset with a new record added for each group (as defined by
`by_vars` parameter). That is, a variable will only be populated in this
new record if it is specified in `by_vars`. For each new record,

- `set_values_to` lists each specified variable and computes its value,

- the variable(s) specified on the LHS of `set_values_to` are set to
  their paired value (RHS). In addition, the start and end date are
  computed as the minimum/maximum dates by `by_vars`.

If the input datasets contains

- both `AxxDTM` and `AxxDT` then all `ASTDTM`,`AENDTM`, `ASTDT`, `AENDT`
  are computed

- only `AxxDTM` then `ASTDTM`,`AENDTM` are computed

- only `AxxDT` then `ASTDT`,`AENDT` are computed.

## Details

For each group (with respect to the variables specified for the
`by_vars` parameter), an observation is added to the output dataset and
the defined values are set to the defined variables

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exist_flag.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
adex <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC, ~VISIT, ~ASTDT, ~AENDT,
  "1015", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
  "1015", "DOSE", 85, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
  "1015", "DOSE", 82, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
  "1015", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
  "1015", "ADJ", NA, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
  "1015", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
  "1017", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
  "1017", "DOSE", 50, NA_character_, "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
  "1017", "DOSE", 65, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02"),
  "1017", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
  "1017", "ADJ", NA, "ADVERSE EVENT", "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
  "1017", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02")
) %>%
  mutate(ASTDTM = ymd_hms(paste(ASTDT, "00:00:00")), AENDTM = ymd_hms(paste(AENDT, "00:00:00")))

# Cumulative dose
adex %>%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(
      PARAMCD = "TDOSE",
      PARCAT1 = "OVERALL",
      AVAL = sum(AVAL, na.rm = TRUE)
    ),
    input_code = "DOSE"
  ) %>%
  select(-ASTDTM, -AENDTM)
#> # A tibble: 14 × 8
#>    USUBJID PARAMCD  AVAL AVALC         VISIT    ASTDT      AENDT      PARCAT1
#>    <chr>   <chr>   <dbl> <chr>         <chr>    <date>     <date>     <chr>  
#>  1 1015    DOSE       80 NA            BASELINE 2014-01-02 2014-01-16 NA     
#>  2 1015    DOSE       85 NA            WEEK 2   2014-01-17 2014-06-18 NA     
#>  3 1015    DOSE       82 NA            WEEK 24  2014-06-19 2014-07-02 NA     
#>  4 1015    ADJ        NA NA            BASELINE 2014-01-02 2014-01-16 NA     
#>  5 1015    ADJ        NA NA            WEEK 2   2014-01-17 2014-06-18 NA     
#>  6 1015    ADJ        NA NA            WEEK 24  2014-06-19 2014-07-02 NA     
#>  7 1017    DOSE       80 NA            BASELINE 2014-01-05 2014-01-19 NA     
#>  8 1017    DOSE       50 NA            WEEK 2   2014-01-20 2014-05-10 NA     
#>  9 1017    DOSE       65 NA            WEEK 24  2014-05-10 2014-07-02 NA     
#> 10 1017    ADJ        NA NA            BASELINE 2014-01-05 2014-01-19 NA     
#> 11 1017    ADJ        NA ADVERSE EVENT WEEK 2   2014-01-20 2014-05-10 NA     
#> 12 1017    ADJ        NA NA            WEEK 24  2014-05-10 2014-07-02 NA     
#> 13 1015    TDOSE     247 NA            NA       2014-01-02 2014-07-02 OVERALL
#> 14 1017    TDOSE     195 NA            NA       2014-01-05 2014-07-02 OVERALL

# average dose in w2-24
adex %>%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    filter_add = VISIT %in% c("WEEK 2", "WEEK 24"),
    set_values_to = exprs(
      PARAMCD = "AVDW224",
      PARCAT1 = "WEEK2-24",
      AVAL = mean(AVAL, na.rm = TRUE)
    ),
    input_code = "DOSE"
  ) %>%
  select(-ASTDTM, -AENDTM)
#> # A tibble: 14 × 8
#>    USUBJID PARAMCD  AVAL AVALC         VISIT    ASTDT      AENDT      PARCAT1 
#>    <chr>   <chr>   <dbl> <chr>         <chr>    <date>     <date>     <chr>   
#>  1 1015    DOSE     80   NA            BASELINE 2014-01-02 2014-01-16 NA      
#>  2 1015    DOSE     85   NA            WEEK 2   2014-01-17 2014-06-18 NA      
#>  3 1015    DOSE     82   NA            WEEK 24  2014-06-19 2014-07-02 NA      
#>  4 1015    ADJ      NA   NA            BASELINE 2014-01-02 2014-01-16 NA      
#>  5 1015    ADJ      NA   NA            WEEK 2   2014-01-17 2014-06-18 NA      
#>  6 1015    ADJ      NA   NA            WEEK 24  2014-06-19 2014-07-02 NA      
#>  7 1017    DOSE     80   NA            BASELINE 2014-01-05 2014-01-19 NA      
#>  8 1017    DOSE     50   NA            WEEK 2   2014-01-20 2014-05-10 NA      
#>  9 1017    DOSE     65   NA            WEEK 24  2014-05-10 2014-07-02 NA      
#> 10 1017    ADJ      NA   NA            BASELINE 2014-01-05 2014-01-19 NA      
#> 11 1017    ADJ      NA   ADVERSE EVENT WEEK 2   2014-01-20 2014-05-10 NA      
#> 12 1017    ADJ      NA   NA            WEEK 24  2014-05-10 2014-07-02 NA      
#> 13 1015    AVDW224  83.5 NA            NA       2014-01-17 2014-07-02 WEEK2-24
#> 14 1017    AVDW224  57.5 NA            NA       2014-01-20 2014-07-02 WEEK2-24

# Any dose adjustment?
adex %>%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(
      PARAMCD = "TADJ",
      PARCAT1 = "OVERALL",
      AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
    ),
    input_code = "ADJ"
  ) %>%
  select(-ASTDTM, -AENDTM)
#> # A tibble: 14 × 8
#>    USUBJID PARAMCD  AVAL AVALC         VISIT    ASTDT      AENDT      PARCAT1
#>    <chr>   <chr>   <dbl> <chr>         <chr>    <date>     <date>     <chr>  
#>  1 1015    DOSE       80 NA            BASELINE 2014-01-02 2014-01-16 NA     
#>  2 1015    DOSE       85 NA            WEEK 2   2014-01-17 2014-06-18 NA     
#>  3 1015    DOSE       82 NA            WEEK 24  2014-06-19 2014-07-02 NA     
#>  4 1015    ADJ        NA NA            BASELINE 2014-01-02 2014-01-16 NA     
#>  5 1015    ADJ        NA NA            WEEK 2   2014-01-17 2014-06-18 NA     
#>  6 1015    ADJ        NA NA            WEEK 24  2014-06-19 2014-07-02 NA     
#>  7 1017    DOSE       80 NA            BASELINE 2014-01-05 2014-01-19 NA     
#>  8 1017    DOSE       50 NA            WEEK 2   2014-01-20 2014-05-10 NA     
#>  9 1017    DOSE       65 NA            WEEK 24  2014-05-10 2014-07-02 NA     
#> 10 1017    ADJ        NA NA            BASELINE 2014-01-05 2014-01-19 NA     
#> 11 1017    ADJ        NA ADVERSE EVENT WEEK 2   2014-01-20 2014-05-10 NA     
#> 12 1017    ADJ        NA NA            WEEK 24  2014-05-10 2014-07-02 NA     
#> 13 1015    TADJ       NA NA            NA       2014-01-02 2014-07-02 OVERALL
#> 14 1017    TADJ       NA Y             NA       2014-01-05 2014-07-02 OVERALL
```
