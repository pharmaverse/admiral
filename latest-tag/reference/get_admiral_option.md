# Get the Value of an Admiral Option

Get the Value of an Admiral Option Which Can Be Modified for Advanced
Users.

## Usage

``` r
get_admiral_option(option)
```

## Arguments

- option:

  A character scalar of commonly used admiral function inputs.

  As of now, support only available for "subject_keys", "signif_digits",
  and "save_memory". See
  [`set_admiral_options()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/set_admiral_options.md)
  for a description of the options.

  Default value

  :   none

## Value

The value of the specified option.

## Details

This function allows flexibility for function inputs that may need to be
repeated multiple times in a script, such as `subject_keys`.

## See also

[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/set_admiral_options.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exist_flag.md),
[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_tte.md)
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_extreme_dtm.md),
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md),
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md)

Other admiral_options:
[`set_admiral_options()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/set_admiral_options.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)
dm <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",   "DM", "01-1302",   61, "YEARS",
  "PILOT01",   "DM", "17-1344",   64, "YEARS"
)

vs <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~VSTESTCD,     ~VISIT,     ~VSTPT, ~VSSTRESN,
  "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE",    "LYING",        76,
  "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE", "STANDING",        87,
  "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2",    "LYING",        71,
  "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2", "STANDING",        79,
  "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE",    "LYING",        88,
  "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE", "STANDING",        86,
  "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2",    "LYING",        84,
  "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2", "STANDING",        82
)

# Merging all dm variables to vs
derive_vars_merged(
  vs,
  dataset_add = select(dm, -DOMAIN),
  by_vars = get_admiral_option("subject_keys")
)
#> # A tibble: 8 × 9
#>   STUDYID DOMAIN USUBJID VSTESTCD VISIT    VSTPT    VSSTRESN   AGE AGEU 
#>   <chr>   <chr>  <chr>   <chr>    <chr>    <chr>       <dbl> <dbl> <chr>
#> 1 PILOT01 VS     01-1302 DIABP    BASELINE LYING          76    61 YEARS
#> 2 PILOT01 VS     01-1302 DIABP    BASELINE STANDING       87    61 YEARS
#> 3 PILOT01 VS     01-1302 DIABP    WEEK 2   LYING          71    61 YEARS
#> 4 PILOT01 VS     01-1302 DIABP    WEEK 2   STANDING       79    61 YEARS
#> 5 PILOT01 VS     17-1344 DIABP    BASELINE LYING          88    64 YEARS
#> 6 PILOT01 VS     17-1344 DIABP    BASELINE STANDING       86    64 YEARS
#> 7 PILOT01 VS     17-1344 DIABP    WEEK 2   LYING          84    64 YEARS
#> 8 PILOT01 VS     17-1344 DIABP    WEEK 2   STANDING       82    64 YEARS
```
