# Derive On-Treatment Flag Variable

Derive on-treatment flag (`ONTRTFL`) in an ADaM dataset with a single
assessment date (e.g `ADT`) or event start and end dates (e.g.
`ASTDT`/`AENDT`).

## Usage

``` r
derive_var_ontrtfl(
  dataset,
  new_var = ONTRTFL,
  start_date,
  end_date = NULL,
  ref_start_date,
  ref_end_date = NULL,
  ref_end_window = 0,
  ignore_time_for_ref_end_date = TRUE,
  filter_pre_timepoint = NULL,
  span_period = FALSE
)
```

## Arguments

- dataset:

  Input dataset

  Required columns are `start_date`, `end_date`, `ref_start_date` and
  `ref_end_date`.

  Default value

  :   none

- new_var:

  On-treatment flag variable name to be created.

  Default value

  :   `ONTRTFL`

- start_date:

  The start date (e.g. `AESDT`) or assessment date (e.g. `ADT`)
  Required; A date or date-time object column is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   none

- end_date:

  The end date of assessment/event (e.g. `AENDT`) A date or date-time
  object column is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Optional; Default is null. If the used and date value is missing on an
  observation, it is assumed the medication is ongoing and `ONTRTFL` is
  set to `"Y"`.

  Default value

  :   `NULL`

- ref_start_date:

  The lower bound of the on-treatment period Required; A date or
  date-time object column is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  Default value

  :   none

- ref_end_date:

  The upper bound of the on-treatment period A date or date-time object
  column is expected.

  Refer to
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_dt.md)
  to impute and derive a date from a date character vector to a date
  object.

  If set to `NULL`, everything after `ref_start_date` will be considered
  on-treatment.

  Default value

  :   `NULL`

- ref_end_window:

  A window to add to the upper bound `ref_end_date` measured in days
  (e.g. 7 if 7 days should be added to the upper bound)

  Default value

  :   `0`

- ignore_time_for_ref_end_date:

  If the argument is set to `TRUE`, the time part is ignored for
  checking if the event occurred more than `ref_end_window` days after
  reference end date.

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `TRUE`

- filter_pre_timepoint:

  An expression to filter observations as not on-treatment when `date` =
  `ref_start_date`. For example, if observations where `VSTPT = PRE`
  should not be considered on-treatment when `date = ref_start_date`,
  `filter_pre_timepoint` should be used to denote when the on-treatment
  flag should be set to null. Optional; default is `NULL`.

  Default value

  :   `NULL`

- span_period:

  A logical scalar. If `TRUE`, events that started prior to the
  `ref_start_date`and are ongoing or end after the `ref_start_date` are
  flagged as `"Y"`. Optional; default is `FALSE`.

  Default value

  :   `FALSE`

## Value

The input dataset with an additional column named `ONTRTFL` with a value
of `"Y"` or `NA`

## Details

On-Treatment is calculated by determining whether the assessment date or
start/stop dates fall between 2 dates. The following logic is used to
assign on-treatment = `"Y"`:

1.  `start_date` is missing and `ref_start_date`is non-missing

2.  No timepoint filter is provided (`filter_pre_timepoint`) and both
    `start_date` and `ref_start_date` are non-missing and `start_date` =
    `ref_start_date`

3.  A timepoint is provided (`filter_pre_timepoint`) and both
    `start_date` and `ref_start_date` are non-missing and
    `start_date = ref_start_date` and the filter provided in
    `filter_pre_timepoint` is not true.

4.  `ref_end_date` is not provided and `ref_start_date < start_date`

5.  `ref_end_date` is provided and `ref_start_date < start_date` \<=
    `ref_end_date + ref_end_window`.

If the `end_date` is provided and the `end_date` \< ref_start_date then
the `ONTRTFL` is set to `NULL`.This would be applicable to cases where
the `start_date` is missing and `ONTRTFL` has been assigned as `"Y"`
above.

If the `span_period` is `TRUE`, this allows the user to assign `ONTRTFL`
as `"Y"` to cases where the record started prior to the `ref_start_date`
and was ongoing or ended after the `ref_start_date`.

Any date imputations needed should be done prior to calling this
function.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_nfrlt.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)

advs <- tribble(
  ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
  "P01",    ymd("2020-02-24"), ymd("2020-01-01"), ymd("2020-03-01"),
  "P02",    ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-03-01"),
  "P03",    ymd("2019-12-31"), ymd("2020-01-01"), ymd("2020-03-01")
)
derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT
)
#> # A tibble: 3 × 5
#>   USUBJID ADT        TRTSDT     TRTEDT     ONTRTFL
#>   <chr>   <date>     <date>     <date>     <chr>  
#> 1 P01     2020-02-24 2020-01-01 2020-03-01 Y      
#> 2 P02     2020-01-01 2020-01-01 2020-03-01 Y      
#> 3 P03     2019-12-31 2020-01-01 2020-03-01 NA     

advs <- tribble(
  ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
  "P01",    ymd("2020-07-01"), ymd("2020-01-01"), ymd("2020-03-01"),
  "P02",    ymd("2020-04-30"), ymd("2020-01-01"), ymd("2020-03-01"),
  "P03",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01")
)
derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  ref_end_window = 60
)
#> # A tibble: 3 × 5
#>   USUBJID ADT        TRTSDT     TRTEDT     ONTRTFL
#>   <chr>   <date>     <date>     <date>     <chr>  
#> 1 P01     2020-07-01 2020-01-01 2020-03-01 NA     
#> 2 P02     2020-04-30 2020-01-01 2020-03-01 Y      
#> 3 P03     2020-03-15 2020-01-01 2020-03-01 Y      

advs <- tribble(
  ~USUBJID, ~ADTM,                      ~TRTSDTM,                   ~TRTEDTM,
  "P01",    ymd_hm("2020-01-02T12:00"), ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
  "P02",    ymd("2020-01-01"),          ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
  "P03",    ymd("2019-12-31"),          ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
) %>%
  mutate(TPT = c(NA, "PRE", NA))
derive_var_ontrtfl(
  advs,
  start_date = ADTM,
  ref_start_date = TRTSDTM,
  ref_end_date = TRTEDTM,
  filter_pre_timepoint = TPT == "PRE"
)
#> # A tibble: 3 × 6
#>   USUBJID ADTM                TRTSDTM             TRTEDTM             TPT  
#>   <chr>   <dttm>              <dttm>              <dttm>              <chr>
#> 1 P01     2020-01-02 12:00:00 2020-01-01 12:00:00 2020-03-01 12:00:00 NA   
#> 2 P02     2020-01-01 00:00:00 2020-01-01 12:00:00 2020-03-01 12:00:00 PRE  
#> 3 P03     2019-12-31 00:00:00 2020-01-01 12:00:00 2020-03-01 12:00:00 NA   
#> # ℹ 1 more variable: ONTRTFL <chr>

advs <- tribble(
  ~USUBJID, ~ASTDT,            ~TRTSDT,           ~TRTEDT,           ~AENDT,
  "P01",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"),
  "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"),
  "P03",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
)
derive_var_ontrtfl(
  advs,
  start_date = ASTDT,
  end_date = AENDT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  ref_end_window = 60,
  span_period = TRUE
)
#> # A tibble: 3 × 6
#>   USUBJID ASTDT      TRTSDT     TRTEDT     AENDT      ONTRTFL
#>   <chr>   <date>     <date>     <date>     <date>     <chr>  
#> 1 P01     2020-03-15 2020-01-01 2020-03-01 2020-12-01 Y      
#> 2 P02     2019-04-30 2020-01-01 2020-03-01 2020-03-15 Y      
#> 3 P03     2019-04-30 2020-01-01 2020-03-01 NA         Y      

advs <- tribble(
  ~USUBJID, ~ASTDT,            ~AP01SDT,          ~AP01EDT,          ~AENDT,
  "P01",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"),
  "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"),
  "P03",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
)
derive_var_ontrtfl(
  advs,
  new_var = ONTR01FL,
  start_date = ASTDT,
  end_date = AENDT,
  ref_start_date = AP01SDT,
  ref_end_date = AP01EDT,
  span_period = TRUE
)
#> # A tibble: 3 × 6
#>   USUBJID ASTDT      AP01SDT    AP01EDT    AENDT      ONTR01FL
#>   <chr>   <date>     <date>     <date>     <date>     <chr>   
#> 1 P01     2020-03-15 2020-01-01 2020-03-01 2020-12-01 NA      
#> 2 P02     2019-04-30 2020-01-01 2020-03-01 2020-03-15 Y       
#> 3 P03     2019-04-30 2020-01-01 2020-03-01 NA         Y       
```
