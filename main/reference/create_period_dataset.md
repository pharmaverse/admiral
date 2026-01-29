# Create a Reference Dataset for Subperiods, Periods, or Phases

The function creates a reference dataset for subperiods, periods, or
phases from the `ADSL` dataset. The reference dataset can be used to
derive subperiod, period, or phase variables like `ASPER`, `ASPRSDT`,
`ASPREDT`, `APERIOD`, `APERSDT`, `APEREDT`, `TRTA`, `APHASEN`, `PHSDTM`,
`PHEDTM`, ... in OCCDS and BDS datasets.

## Usage

``` r
create_period_dataset(
  dataset,
  new_vars,
  subject_keys = get_admiral_option("subject_keys")
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `new_vars` and `subject_keys` arguments
  are expected to be in the dataset. For each element of `new_vars` at
  least one variable of the form of the right hand side value must be
  available in the dataset.

  Default value

  :   none

- new_vars:

  New variables

  A named list of variables like
  `exprs(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)` is expected.
  The left hand side of the elements defines a variable of the output
  dataset, the right hand side defines the source variables from the
  ADSL dataset in CDISC notation.

  If the lower case letter "w" is used it refers to a phase variable, if
  the lower case letters "xx" are used it refers to a period variable,
  and if both "xx" and "w" are used it refers to a subperiod variable.

  Only one type must be used, e.g., all right hand side values must
  refer to period variables. It is not allowed to mix for example period
  and subperiod variables. If period *and* subperiod variables are
  required, separate reference datasets must be created.

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of expressions where the expressions are symbols as returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  is expected.

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

A period reference dataset (see "Details" section)

## Details

For each subject and each subperiod/period/phase where at least one of
the source variable is not `NA` an observation is added to the output
dataset.

Depending on the type of the source variable (subperiod, period, or
phase) the variable `ASPER`, `APERIOD`, or `APHASEN` is added and set to
the number of the subperiod, period, or phase.

The variables specified for `new_vars` (left hand side) are added to the
output dataset and set to the value of the source variable (right hand
side).

## See also

[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_period.md)

Creating auxiliary datasets:
[`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/main/reference/consolidate_metadata.md),
[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md),
[`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

# Create reference dataset for periods
adsl <- tribble(
  ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,     ~TRT01A, ~TRT02A,
  "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "A",     "B",
  "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01", "B",     "A",
) %>%
  mutate(
    across(matches("AP\\d\\d[ES]DT"), ymd)
  ) %>%
  mutate(
    STUDYID = "xyz"
  )

create_period_dataset(
  adsl,
  new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA)
)
#> # A tibble: 4 × 6
#>   STUDYID USUBJID APERIOD APERSDT    APEREDT    TRTA 
#>   <chr>   <chr>     <int> <date>     <date>     <chr>
#> 1 xyz     1             1 2021-01-04 2021-02-06 A    
#> 2 xyz     1             2 2021-02-07 2021-03-07 B    
#> 3 xyz     2             1 2021-02-02 2021-03-02 B    
#> 4 xyz     2             2 2021-03-03 2021-04-01 A    

# Create reference dataset for phases
adsl <- tribble(
  ~USUBJID, ~PH1SDT,      ~PH1EDT,      ~PH2SDT,      ~PH2EDT,      ~APHASE1,    ~APHASE2,
  "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "TREATMENT", "FUP",
  "2",      "2021-02-02", "2021-03-02", NA,           NA,           "TREATMENT", NA
) %>%
  mutate(
    across(matches("PH\\d[ES]DT"), ymd)
  ) %>%
  mutate(
    STUDYID = "xyz"
  )

create_period_dataset(
  adsl,
  new_vars = exprs(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
)
#> # A tibble: 3 × 6
#>   STUDYID USUBJID APHASEN PHSDT      PHEDT      APHASE   
#>   <chr>   <chr>     <int> <date>     <date>     <chr>    
#> 1 xyz     1             1 2021-01-04 2021-02-06 TREATMENT
#> 2 xyz     1             2 2021-02-07 2021-03-07 FUP      
#> 3 xyz     2             1 2021-02-02 2021-03-02 TREATMENT

# Create reference datasets for subperiods
adsl <- tribble(
  ~USUBJID, ~P01S1SDT,    ~P01S1EDT,    ~P01S2SDT,    ~P01S2EDT,    ~P02S1SDT,    ~P02S1EDT,
  "1",      "2021-01-04", "2021-01-19", "2021-01-20", "2021-02-06", "2021-02-07", "2021-03-07",
  "2",      "2021-02-02", "2021-03-02", NA,           NA,           "2021-03-03", "2021-04-01"
) %>%
  mutate(
    across(matches("P\\d\\dS\\d[ES]DT"), ymd)
  ) %>%
  mutate(
    STUDYID = "xyz"
  )

create_period_dataset(
  adsl,
  new_vars = exprs(ASPRSDT = PxxSwSDT, ASPREDT = PxxSwEDT)
)
#> # A tibble: 5 × 6
#>   STUDYID USUBJID APERIOD ASPER ASPRSDT    ASPREDT   
#>   <chr>   <chr>     <int> <int> <date>     <date>    
#> 1 xyz     1             1     1 2021-01-04 2021-01-19
#> 2 xyz     1             1     2 2021-01-20 2021-02-06
#> 3 xyz     1             2     1 2021-02-07 2021-03-07
#> 4 xyz     2             1     1 2021-02-02 2021-03-02
#> 5 xyz     2             2     1 2021-03-03 2021-04-01
```
