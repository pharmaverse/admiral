# Derive First or Last Date from Multiple Sources

**\[deprecated\]** The `derive_var_extreme_dt()` function has been
deprecated in favor of
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_extreme_event.md).

Add the first or last date from multiple sources to the dataset, e.g.,
the last known alive date (`LSTALVDT`).

**Note:** This is a wrapper function for the function
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md).

## Usage

``` r
derive_var_extreme_dt(
  dataset,
  new_var,
  ...,
  source_datasets,
  mode,
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

- new_var:

  Name of variable to create

  Default value

  :   none

- ...:

  Source(s) of dates. One or more
  [`date_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/date_source.md)
  objects are expected.

  Default value

  :   none

- source_datasets:

  A named `list` containing datasets in which to search for the first or
  last date

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, the first date for each subject is
  selected. If `"last"` is specified, the last date for each subject is
  selected.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- subject_keys:

  Variables to uniquely identify a subject

  A list of expressions where the expressions are symbols as returned by
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
  is expected.

  Default value

  :   `get_admiral_option("subject_keys")`

## Value

The input dataset with the new variable added.

## Details

The following steps are performed to create the output dataset:

1.  For each source dataset the observations as specified by the
    `filter` element are selected and observations where `date` is `NA`
    are removed. Then for each patient the first or last observation
    (with respect to `date` and `mode`) is selected.

2.  The new variable is set to the variable or expression specified by
    the `date` element.

3.  The variables specified by the `set_values_to` element are added.

4.  The selected observations of all source datasets are combined into a
    single dataset.

5.  For each patient the first or last observation (with respect to the
    new variable and `mode`) from the single dataset is selected and the
    new variable is merged to the input dataset.

6.  The time part is removed from the new variable.

## See also

[`date_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/date_source.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md),
[`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_merged.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_extreme_record.md),
[`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_dthcaus.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/get_summary_records.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)
ae <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AESEQ,     ~AESTDTC,     ~AEENDTC,
  "PILOT01",    "AE", "01-1130",      5, "2014-05-09", "2014-05-09",
  "PILOT01",    "AE", "01-1130",      6, "2014-05-22",           NA,
  "PILOT01",    "AE", "01-1130",      4, "2014-05-09", "2014-05-09",
  "PILOT01",    "AE", "01-1130",      8, "2014-05-22",           NA,
  "PILOT01",    "AE", "01-1130",      7, "2014-05-22",           NA,
  "PILOT01",    "AE", "01-1130",      2, "2014-03-09", "2014-03-09",
  "PILOT01",    "AE", "01-1130",      1, "2014-03-09", "2014-03-16",
  "PILOT01",    "AE", "01-1130",      3, "2014-03-09", "2014-03-16",
  "PILOT01",    "AE", "01-1133",      1, "2012-12-27",           NA,
  "PILOT01",    "AE", "01-1133",      3, "2012-12-27",           NA,
  "PILOT01",    "AE", "01-1133",      2, "2012-12-27",           NA,
  "PILOT01",    "AE", "01-1133",      4, "2012-12-27",           NA,
  "PILOT01",    "AE", "01-1211",      5, "2012-11-29",           NA,
  "PILOT01",    "AE", "01-1211",      1, "2012-11-16",           NA,
  "PILOT01",    "AE", "01-1211",      7, "2013-01-11",           NA,
  "PILOT01",    "AE", "01-1211",      8, "2013-01-11",           NA,
  "PILOT01",    "AE", "01-1211",      4, "2012-11-22",           NA,
  "PILOT01",    "AE", "01-1211",      2, "2012-11-21", "2012-11-21",
  "PILOT01",    "AE", "01-1211",      3, "2012-11-21",           NA,
  "PILOT01",    "AE", "01-1211",      6, "2012-12-09",           NA,
  "PILOT01",    "AE", "01-1211",      9, "2013-01-14", "2013-01-14",
  "PILOT01",    "AE", "09-1081",      2, "2014-05-01",           NA,
  "PILOT01",    "AE", "09-1081",      1, "2014-04-07",           NA,
  "PILOT01",    "AE", "09-1088",      1, "2014-05-08",           NA,
  "PILOT01",    "AE", "09-1088",      2, "2014-08-02",           NA
)

adsl <- tribble(
  ~STUDYID,   ~USUBJID,              ~TRTEDTM,      ~TRTEDT,
  "PILOT01", "01-1130", "2014-08-16 23:59:59", "2014-08-16",
  "PILOT01", "01-1133", "2013-04-28 23:59:59", "2013-04-28",
  "PILOT01", "01-1211", "2013-01-12 23:59:59", "2013-01-12",
  "PILOT01", "09-1081", "2014-04-27 23:59:59", "2014-04-27",
  "PILOT01", "09-1088", "2014-10-09 23:59:59", "2014-10-09"
) %>%
  mutate(
    across(TRTEDTM:TRTEDT, as.Date)
  )


lb <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
  "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
  "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
  "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
  "PILOT01",    "LB", "01-1133",    304, "2013-04-29T10:13",
  "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
  "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
  "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
  "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
  "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
  "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
)

dm <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1130",   84, "YEARS",
  "PILOT01",    "DM", "01-1133",   81, "YEARS",
  "PILOT01",    "DM", "01-1211",   76, "YEARS",
  "PILOT01",    "DM", "09-1081",   86, "YEARS",
  "PILOT01",    "DM", "09-1088",   69, "YEARS"
)

ae_start <- date_source(
  dataset_name = "ae",
  date = convert_dtc_to_dt(AESTDTC, highest_imputation = "M")
)
ae_end <- date_source(
  dataset_name = "ae",
  date = convert_dtc_to_dt(AEENDTC, highest_imputation = "M")
)

ae_ext <- ae %>%
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AEST",
    highest_imputation = "M"
  ) %>%
  derive_vars_dt(
    dtc = AEENDTC,
    new_vars_prefix = "AEEN",
    highest_imputation = "M"
  )

lb_date <- date_source(
  dataset_name = "lb",
  date = convert_dtc_to_dt(LBDTC)
)

lb_ext <- derive_vars_dt(
  lb,
  dtc = LBDTC,
  new_vars_prefix = "LB"
)

adsl_date <- date_source(dataset_name = "adsl", date = TRTEDT)

dm %>%
  derive_var_extreme_dt(
    new_var = LSTALVDT,
    ae_start, ae_end, lb_date, adsl_date,
    source_datasets = list(
      adsl = adsl,
      ae = ae_ext,
      lb = lb_ext
    ),
    mode = "last"
  ) %>%
  select(USUBJID, LSTALVDT)
#> Warning: `derive_var_extreme_dt()` was deprecated in admiral 1.2.0.
#> ℹ Please use `derive_vars_extreme_event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
#> Warning: `derive_var_extreme_dtm()` was deprecated in admiral 1.2.0.
#> ℹ Please use `derive_vars_extreme_event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
#> ℹ The deprecated feature was likely used in the admiral package.
#>   Please report the issue at <https://github.com/pharmaverse/admiral/issues>.
#> # A tibble: 5 × 2
#>   USUBJID LSTALVDT  
#>   <chr>   <date>    
#> 1 01-1130 2014-08-16
#> 2 01-1133 2013-04-29
#> 3 01-1211 2013-01-14
#> 4 09-1081 2014-05-10
#> 5 09-1088 2014-10-09

# derive last alive date and traceability variables
ae_start <- date_source(
  dataset_name = "ae",
  date = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
  set_values_to = exprs(
    LALVDOM = "AE",
    LALVSEQ = AESEQ,
    LALVVAR = "AESTDTC"
  )
)

ae_end <- date_source(
  dataset_name = "ae",
  date = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
  set_values_to = exprs(
    LALVDOM = "AE",
    LALVSEQ = AESEQ,
    LALVVAR = "AEENDTC"
  )
)

lb_date <- date_source(
  dataset_name = "lb",
  date = convert_dtc_to_dt(LBDTC),
  set_values_to = exprs(
    LALVDOM = "LB",
    LALVSEQ = LBSEQ,
    LALVVAR = "LBDTC"
  )
)

adsl_date <- date_source(
  dataset_name = "adsl",
  date = TRTEDT,
  set_values_to = exprs(
    LALVDOM = "ADSL",
    LALVSEQ = NA_integer_,
    LALVVAR = "TRTEDT"
  )
)

dm %>%
  derive_var_extreme_dt(
    new_var = LSTALVDT,
    ae_start, ae_end, lb_date, adsl_date,
    source_datasets = list(
      adsl = adsl,
      ae = ae_ext,
      lb = lb_ext
    ),
    mode = "last"
  ) %>%
  select(USUBJID, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR)
#> # A tibble: 5 × 5
#>   USUBJID LSTALVDT   LALVDOM LALVSEQ LALVVAR
#>   <chr>   <date>     <chr>     <dbl> <chr>  
#> 1 01-1130 2014-08-16 ADSL         NA TRTEDT 
#> 2 01-1133 2013-04-29 LB          304 LBDTC  
#> 3 01-1211 2013-01-14 AE            9 AEENDTC
#> 4 09-1081 2014-05-10 LB          219 LBDTC  
#> 5 09-1088 2014-10-09 ADSL         NA TRTEDT 
```
