# Set the Value of admiral Options

Set the values of admiral options that can be modified for advanced
users.

## Usage

``` r
set_admiral_options(subject_keys, signif_digits, save_memory)
```

## Arguments

- subject_keys:

  Variables to uniquely identify a subject, defaults to
  `exprs(STUDYID, USUBJID)`. This option is used as default value for
  the `subject_keys` argument in all admiral functions.

  Default value

  :   none

- signif_digits:

  Holds number of significant digits when comparing to numeric
  variables, defaults to `15`. This option is used as default value for
  the `signif_dig` argument in admiral functions
  [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  and
  [`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md).

  Default value

  :   none

- save_memory:

  If set to `TRUE`, an alternative algorithm is used in the functions
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
  and
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  which requires less memory but more run-time.

  Default value

  :   none

## Value

No return value, called for side effects.

## Details

Modify an admiral option, e.g `subject_keys`, such that it automatically
affects downstream function inputs where
[`get_admiral_option()`](https:/pharmaverse.github.io/admiral/main/reference/get_admiral_option.md)
is called such as
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md).

## See also

[`get_admiral_option()`](https:/pharmaverse.github.io/admiral/main/reference/get_admiral_option.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_period.md),
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_period_dataset.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md)

Other admiral_options:
[`get_admiral_option()`](https:/pharmaverse.github.io/admiral/main/reference/get_admiral_option.md)

## Examples

``` r
library(lubridate)
library(dplyr, warn.conflicts = FALSE)
library(tibble)
set_admiral_options(subject_keys = exprs(STUDYID, USUBJID2))

# Derive a new parameter for measurable disease at baseline
adsl <- tribble(
  ~USUBJID2,
  "1",
  "2",
  "3"
) %>%
  mutate(STUDYID = "XX1234")

tu <- tribble(
  ~USUBJID2,      ~VISIT,    ~TUSTRESC,
  "1",       "SCREENING",     "TARGET",
  "1",          "WEEK 1",     "TARGET",
  "1",          "WEEK 5",     "TARGET",
  "1",          "WEEK 9", "NON-TARGET",
  "2",       "SCREENING", "NON-TARGET",
  "2",       "SCREENING", "NON-TARGET"
) %>%
  mutate(
    STUDYID = "XX1234",
    TUTESTCD = "TUMIDENT"
  )

derive_param_exist_flag(
  dataset_ref = adsl,
  dataset_add = tu,
  filter_add = TUTESTCD == "TUMIDENT" & VISIT == "SCREENING",
  condition = TUSTRESC == "TARGET",
  false_value = "N",
  missing_value = "N",
  set_values_to = exprs(
    PARAMCD = "MDIS",
    PARAM = "Measurable Disease at Baseline"
  )
)
#> # A tibble: 3 × 5
#>   STUDYID USUBJID2 AVALC PARAMCD PARAM                         
#>   <chr>   <chr>    <chr> <chr>   <chr>                         
#> 1 XX1234  1        Y     MDIS    Measurable Disease at Baseline
#> 2 XX1234  2        N     MDIS    Measurable Disease at Baseline
#> 3 XX1234  3        N     MDIS    Measurable Disease at Baseline

set_admiral_options(signif_digits = 14)

# Derive ANRIND for ADVS
advs <- tribble(
  ~PARAMCD, ~AVAL, ~ANRLO, ~ANRHI,
  "DIABP",     59,     60,     80,
  "SYSBP",    120,     90,    130,
  "RESP",      21,      8,     20,
)

derive_var_anrind(advs)
#> # A tibble: 3 × 5
#>   PARAMCD  AVAL ANRLO ANRHI ANRIND
#>   <chr>   <dbl> <dbl> <dbl> <chr> 
#> 1 DIABP      59    60    80 LOW   
#> 2 SYSBP     120    90   130 NORMAL
#> 3 RESP       21     8    20 HIGH  
```
