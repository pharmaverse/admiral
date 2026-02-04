# Select the First or Last Date from Several Sources

Select for each subject the first or last observation with respect to a
date from a list of sources.

## Usage

``` r
filter_date_sources(
  sources,
  source_datasets,
  by_vars,
  create_datetime = FALSE,
  subject_keys,
  mode,
  check_type = "none"
)
```

## Arguments

- sources:

  Sources

  A list of
  [`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)
  objects is expected.

  Default value

  :   none

- source_datasets:

  Source datasets

  A named list of datasets is expected. The `dataset_name` field of
  [`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)
  refers to the dataset provided in the list.

  Default value

  :   none

- by_vars:

  By variables

  If the parameter is specified, for each by group the observations are
  selected separately.

  Default value

  :   none

- create_datetime:

  Create datetime variable?

  If set to `TRUE`, variables `ADTM` is created. Otherwise, variables
  `ADT` is created.

  Default value

  :   `FALSE`

- subject_keys:

  Variables to uniquely identify a subject

  A list of symbols created using
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  is expected.

  Default value

  :   none

- mode:

  Selection mode (first or last)

  If `"first"` is specified, for each subject the first observation with
  respect to the date is included in the output dataset. If `"last"` is
  specified, the last observation is included in the output dataset.

  Permitted values

  :   `"first"`, `"last"`

  Default value

  :   none

- check_type:

  Check uniqueness

  If `"warning"`, `"message"`, or `"error"` is specified, the specified
  message is issued if the observations of the source datasets are not
  unique with respect to the by variables and the date and order
  specified in the
  [`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)
  objects.

  Permitted values

  :   `"none"`, `"warning"`, `"error"`, `"message"`

  Default value

  :   `"none"`

## Value

A dataset with one observation per subject as described in the "Details"
section.

## Details

The following steps are performed to create the output dataset:

1.  For each source dataset the observations as specified by the
    `filter` element are selected. Then for each subject the first or
    last observation (with respect to `date`) is selected.

2.  The `ADT` variable is set to the variable specified by the `date`
    element. If the date variable is a datetime variable, only the
    datepart is copied. If the source variable is a character variable,
    it is converted to a date. If the date is incomplete, it is imputed
    as the first possible date.

3.  The `CNSR` is added and set to the value of the `censor` element.

4.  The selected observations of all source datasets are combined into a
    single dataset.

5.  For each subject the first or last observation (with respect to the
    `ADT` variable) from the single dataset is selected.

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

adsl <- tribble(
  ~USUBJID, ~TRTSDT,           ~EOSDT,
  "01",     ymd("2020-12-06"), ymd("2021-03-06"),
  "02",     ymd("2021-01-16"), ymd("2021-02-03")
) %>%
  mutate(STUDYID = "AB42")

ae <- tribble(
  ~USUBJID, ~AESTDTC,     ~AESEQ, ~AEDECOD,
  "01",     "2021-01-03", 1,      "Flu",
  "01",     "2021-03-04", 2,      "Cough",
  "01",     "2021-01-01", 3,      "Flu"
) %>%
  mutate(
    STUDYID = "AB42",
    AESTDT = ymd(AESTDTC)
  )

ttae <- event_source(
  dataset_name = "ae",
  date = AESTDT,
  set_values_to = exprs(
    EVNTDESC = "AE",
    SRCDOM = "AE",
    SRCVAR = "AESTDTC",
    SRCSEQ = AESEQ
  )
)

admiral:::filter_date_sources(
  sources = list(ttae),
  source_datasets = list(adsl = adsl, ae = ae),
  by_vars = exprs(AEDECOD),
  create_datetime = FALSE,
  subject_keys = get_admiral_option("subject_keys"),
  mode = "first",
  check_type = "none"
)
#> # A tibble: 2 × 9
#>   USUBJID AEDECOD STUDYID EVNTDESC SRCDOM SRCVAR  SRCSEQ  CNSR ADT       
#>   <chr>   <chr>   <chr>   <chr>    <chr>  <chr>    <dbl> <int> <date>    
#> 1 01      Cough   AB42    AE       AE     AESTDTC      2     0 2021-03-04
#> 2 01      Flu     AB42    AE       AE     AESTDTC      3     0 2021-01-01
```
