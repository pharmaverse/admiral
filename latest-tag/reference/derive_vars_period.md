# Add Subperiod, Period, or Phase Variables to ADSL

The function adds subperiod, period, or phase variables like `P01S1SDT`,
`P01S2SDT`, `AP01SDTM`, `AP02SDTM`, `TRT01A`, `TRT02A`, `PH1SDT`,
`PH2SDT`, ... to the input dataset. The values of the variables are
defined by a period reference dataset which has one observations per
patient and subperiod, period, or phase.

## Usage

``` r
derive_vars_period(
  dataset,
  dataset_ref,
  new_vars,
  subject_keys = get_admiral_option("subject_keys")
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `subject_keys` argument are expected to
  be in the dataset.

  Default value

  :   none

- dataset_ref:

  Period reference dataset

  The variables specified by `new_vars` and `subject_keys` are expected.

  If subperiod variables are requested, `APERIOD` and `ASPER` are
  expected. If period variables are requested. `APERIOD` is expected. If
  phase variables are requested, `APHASEN` is expected.

  Default value

  :   none

- new_vars:

  New variables

  A named list of variables like
  `exprs(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)` is expected.
  The left hand side of the elements defines a set of variables (in
  CDISC notation) to be added to the output dataset. The right hand side
  defines the source variable from the period reference dataset.

  If the lower case letter "w" is used it refers to a phase variable, if
  the lower case letters "xx" are used it refers to a period variable,
  and if both "xx" and "w" are used it refers to a subperiod variable.

  Only one type must be used, e.g., all left hand side values must refer
  to period variables. It is not allowed to mix for example period and
  subperiod variables. If period *and* subperiod variables are required,
  separate calls must be used.

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of expressions where the expressions are symbols as returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/reexport-exprs.md)
  is expected.

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

The input dataset with subperiod/period/phase variables added (see
"Details" section)

## Details

For each subperiod/period/phase in the period reference dataset and each
element in `new_vars` a variable (LHS value of `new_vars`) is added to
the output dataset and set to the value of the source variable (RHS
value of `new_vars`.

## See also

[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md)

ADSL Functions that returns variable appended to dataset:
[`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_age_years.md),
[`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_aage.md),
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

# Add period variables to ADSL
period_ref <- tribble(
  ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
  "1",             1, "2021-01-04", "2021-02-06",
  "1",             2, "2021-02-07", "2021-03-07",
  "2",             1, "2021-02-02", "2021-03-02",
  "2",             2, "2021-03-03", "2021-04-01"
) %>%
  mutate(
    STUDYID = "xyz",
    APERIOD = as.integer(APERIOD),
    across(matches("APER[ES]DT"), ymd)
  )

derive_vars_period(
  adsl,
  dataset_ref = period_ref,
  new_vars = exprs(APxxSDT = APERSDT, APxxEDT = APEREDT)
) %>%
  select(STUDYID, USUBJID, AP01SDT, AP01EDT, AP02SDT, AP02EDT)
#> # A tibble: 2 × 6
#>   STUDYID USUBJID AP01SDT    AP01EDT    AP02SDT    AP02EDT   
#>   <chr>   <chr>   <date>     <date>     <date>     <date>    
#> 1 xyz     1       2021-01-04 2021-02-06 2021-02-07 2021-03-07
#> 2 xyz     2       2021-02-02 2021-03-02 2021-03-03 2021-04-01

# Add phase variables to ADSL
phase_ref <- tribble(
  ~USUBJID, ~APHASEN, ~PHSDT,       ~PHEDT,       ~APHASE,
  "1",             1, "2021-01-04", "2021-02-06", "TREATMENT",
  "1",             2, "2021-02-07", "2021-03-07", "FUP",
  "2",             1, "2021-02-02", "2021-03-02", "TREATMENT"
) %>%
  mutate(
    STUDYID = "xyz",
    APHASEN = as.integer(APHASEN),
    across(matches("PH[ES]DT"), ymd)
  )

derive_vars_period(
  adsl,
  dataset_ref = phase_ref,
  new_vars = exprs(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)
) %>%
  select(STUDYID, USUBJID, PH1SDT, PH1EDT, PH2SDT, PH2EDT, APHASE1, APHASE2)
#> # A tibble: 2 × 8
#>   STUDYID USUBJID PH1SDT     PH1EDT     PH2SDT     PH2EDT     APHASE1   APHASE2
#>   <chr>   <chr>   <date>     <date>     <date>     <date>     <chr>     <chr>  
#> 1 xyz     1       2021-01-04 2021-02-06 2021-02-07 2021-03-07 TREATMENT FUP    
#> 2 xyz     2       2021-02-02 2021-03-02 NA         NA         TREATMENT NA     

# Add subperiod variables to ADSL
subperiod_ref <- tribble(
  ~USUBJID, ~APERIOD, ~ASPER, ~ASPRSDT,     ~ASPREDT,
  "1",             1,      1, "2021-01-04", "2021-01-19",
  "1",             1,      2, "2021-01-20", "2021-02-06",
  "1",             2,      1, "2021-02-07", "2021-03-07",
  "2",             1,      1, "2021-02-02", "2021-03-02",
  "2",             2,      1, "2021-03-03", "2021-04-01"
) %>%
  mutate(
    STUDYID = "xyz",
    APERIOD = as.integer(APERIOD),
    ASPER = as.integer(ASPER),
    across(matches("ASPR[ES]DT"), ymd)
  )

derive_vars_period(
  adsl,
  dataset_ref = subperiod_ref,
  new_vars = exprs(PxxSwSDT = ASPRSDT, PxxSwEDT = ASPREDT)
) %>%
  select(STUDYID, USUBJID, P01S1SDT, P01S1EDT, P01S2SDT, P01S2EDT, P02S1SDT, P02S1EDT)
#> # A tibble: 2 × 8
#>   STUDYID USUBJID P01S1SDT   P01S1EDT   P01S2SDT   P01S2EDT   P02S1SDT  
#>   <chr>   <chr>   <date>     <date>     <date>     <date>     <date>    
#> 1 xyz     1       2021-01-04 2021-01-19 2021-01-20 2021-02-06 2021-02-07
#> 2 xyz     2       2021-02-02 2021-03-02 NA         NA         2021-03-03
#> # ℹ 1 more variable: P02S1EDT <date>
```
