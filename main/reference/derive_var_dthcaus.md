# Derive Death Cause

**\[deprecated\]** The `derive_var_dthcaus()` function has been
deprecated in favor of
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md).

Derive death cause (`DTHCAUS`) and add traceability variables if
required.

## Usage

``` r
derive_var_dthcaus(
  dataset,
  ...,
  source_datasets,
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

- ...:

  Objects of class "dthcaus_source" created by
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md).

  Default value

  :   none

- source_datasets:

  A named `list` containing datasets in which to search for the death
  cause

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

The input dataset with `DTHCAUS` variable added.

## Details

This function derives `DTHCAUS` along with the user-defined traceability
variables, if required. If a subject has death info from multiple
sources, the one from the source with the earliest death date will be
used. If dates are equivalent, the first source will be kept, so the
user should provide the inputs in the preferred order.

## See also

[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)

Other deprecated:
[`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md),
[`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
[`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md),
[`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md),
[`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),
[`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md),
[`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md),
[`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

adsl <- tribble(
  ~STUDYID,  ~USUBJID,
  "STUDY01", "PAT01",
  "STUDY01", "PAT02",
  "STUDY01", "PAT03"
)
ae <- tribble(
  ~STUDYID,  ~USUBJID, ~AESEQ, ~AEDECOD,       ~AEOUT,  ~AEDTHDTC,
  "STUDY01", "PAT01",  12,     "SUDDEN DEATH", "FATAL", "2021-04-04"
)

ds <- tribble(
  ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
  "STUDY01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
  "STUDY01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
  "STUDY01", "PAT02", 3, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
  "STUDY01", "PAT03", 1, "DEATH", "POST STUDY REPORTING OF DEATH", "2022-03-03"
)

# Derive `DTHCAUS` only - for on-study deaths only
src_ae <- dthcaus_source(
  dataset_name = "ae",
  filter = AEOUT == "FATAL",
  date = convert_dtc_to_dt(AEDTHDTC),
  mode = "first",
  dthcaus = AEDECOD
)
#> Warning: `dthcaus_source()` was deprecated in admiral 1.2.0.
#> ℹ Please use `event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation

src_ds <- dthcaus_source(
  dataset_name = "ds",
  filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
  date = convert_dtc_to_dt(DSSTDTC),
  mode = "first",
  dthcaus = DSTERM
)

derive_var_dthcaus(adsl, src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))
#> Warning: `derive_var_dthcaus()` was deprecated in admiral 1.2.0.
#> ℹ Please use `derive_vars_extreme_event()` instead.
#> ✖ This message will turn into an error at the beginning of 2027.
#> ℹ See admiral's deprecation guidance:
#>   https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
#> # A tibble: 3 × 3
#>   STUDYID USUBJID DTHCAUS                            
#>   <chr>   <chr>   <chr>                              
#> 1 STUDY01 PAT01   SUDDEN DEATH                       
#> 2 STUDY01 PAT02   DEATH DUE TO PROGRESSION OF DISEASE
#> 3 STUDY01 PAT03   NA                                 

# Derive `DTHCAUS` and add traceability variables - for on-study deaths only
src_ae <- dthcaus_source(
  dataset_name = "ae",
  filter = AEOUT == "FATAL",
  date = convert_dtc_to_dt(AEDTHDTC),
  mode = "first",
  dthcaus = AEDECOD,
  set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
)

src_ds <- dthcaus_source(
  dataset_name = "ds",
  filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
  date = convert_dtc_to_dt(DSSTDTC),
  mode = "first",
  dthcaus = DSTERM,
  set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
)

derive_var_dthcaus(adsl, src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))
#> # A tibble: 3 × 5
#>   STUDYID USUBJID DTHCAUS                             DTHDOM DTHSEQ
#>   <chr>   <chr>   <chr>                               <chr>   <dbl>
#> 1 STUDY01 PAT01   SUDDEN DEATH                        AE         12
#> 2 STUDY01 PAT02   DEATH DUE TO PROGRESSION OF DISEASE DS          3
#> 3 STUDY01 PAT03   NA                                  NA         NA

# Derive `DTHCAUS` as above - now including post-study deaths with different `DTHCAUS` value
src_ae <- dthcaus_source(
  dataset_name = "ae",
  filter = AEOUT == "FATAL",
  date = convert_dtc_to_dt(AEDTHDTC),
  mode = "first",
  dthcaus = AEDECOD,
  set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
)

ds <- mutate(
  ds,
  DSSTDT = convert_dtc_to_dt(DSSTDTC)
)

src_ds <- dthcaus_source(
  dataset_name = "ds",
  filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
  date = DSSTDT,
  mode = "first",
  dthcaus = DSTERM,
  set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
)

src_ds_post <- dthcaus_source(
  dataset_name = "ds",
  filter = DSDECOD == "DEATH" & DSTERM == "POST STUDY REPORTING OF DEATH",
  date = DSSTDT,
  mode = "first",
  dthcaus = "POST STUDY: UNKNOWN CAUSE",
  set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
)

derive_var_dthcaus(
  adsl,
  src_ae, src_ds, src_ds_post,
  source_datasets = list(ae = ae, ds = ds)
)
#> # A tibble: 3 × 5
#>   STUDYID USUBJID DTHCAUS                             DTHDOM DTHSEQ
#>   <chr>   <chr>   <chr>                               <chr>   <dbl>
#> 1 STUDY01 PAT01   SUDDEN DEATH                        AE         12
#> 2 STUDY01 PAT02   DEATH DUE TO PROGRESSION OF DISEASE DS          3
#> 3 STUDY01 PAT03   POST STUDY: UNKNOWN CAUSE           DS          1
```
