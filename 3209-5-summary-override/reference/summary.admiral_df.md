# Summarize an `admiral_df` Dataset

Reports what a dataset contains and whether it has one record per key
variable combination. This is useful after a derivation, e.g. to check
that
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/derive_param_computed.md)
added the records that were expected.

## Usage

``` r
# S3 method for class 'admiral_df'
summary(object, adsl = NULL, ...)
```

## Arguments

- object:

  An `admiral_df` object

  Default value

  :   none

- adsl:

  The subject-level dataset to compare `object` with

  ADSL holds every subject of the study, so it shows what `object`
  cannot show by itself: the subjects without any records, and the
  subjects in `object` which should not be there. If it is supplied, the
  summary also reports how many of the ADSL subjects have records,
  overall, within the safety population (`SAFFL`), and per treatment
  arm; checks that every subject of `object` is in ADSL; and, for an
  `OCCDS` dataset, reports the subjects with a treatment-emergent record
  (`TRTEMFL`) out of the safety population, or out of ADSL if there is
  no `SAFFL`.

  Subjects are matched by the variables they are counted by (see
  Details), which must be in `adsl`, and `adsl` must have one record per
  subject. `adsl` is ignored, with a warning, for a dataset of class
  `"ADSL"` and for a dataset without any of these variables.

  Permitted values

  :   a data frame or `NULL`

  Default value

  :   `NULL`

- ...:

  Not used

  Default value

  :   none

## Value

An object of class `summary_admiral_df`, a named list which is printed
by
[`print.summary_admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/print.summary_admiral_df.md).

## Details

The summary always shows the dataset name and ADaM class, the number of
subjects, observations and variables, and the record structure. It then
adds the information for the dataset class (see
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_class.md)):

- `ADSL`: subjects per treatment arm, the population flag counts, the
  subject flow (randomized, treated, completed), the treatment duration,
  the end of study status, the deaths and their causes, and the age and
  sex of the subjects.

- `BDS`/`MDBDS`: records, subjects, visits, and the distribution of
  `AVAL` per parameter, the analysis visits, the derived records per
  `DTYPE`, and a check for more than one baseline record (`ABLFL`) per
  subject and parameter.

- `OCCDS`/`MDOCCDS`: the number of distinct terms, the
  treatment-emergent, and the serious vs total records, the records per
  severity or toxicity grade, and a check for more than one `"Y"` per
  group in the occurrence flags (`AOCCFL`, `AOCCSFL`, `AOCCPFL`, ...).

- `TTE`/`MDTTE`: events and censored per parameter.

If `adsl` is supplied, the summary also compares the dataset with it,
see the `adsl` argument.

The summary reports facts only, apart from the record structure check
and the checks for the dataset class, which report problems without
signaling them. Each item is only reported when the variables it needs
are in the dataset, so a dataset which is still being derived gives a
shorter summary rather than an error.

The record structure is checked against the keys stored by
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).
Keys are never inferred, so if none were stored the structure is not
checked.

Subjects are counted by `USUBJID` if it is present, as the other subject
keys (e.g. `STUDYID`) may not be populated on records added by a
derivation. Otherwise the subject keys
(`get_admiral_option("subject_keys")`) are used if all of them are in
the dataset, as only together do they identify a subject. A record with
any of these missing (`NA` or `""`) belongs to no subject, so it is not
counted as one.

The records of a medical device dataset (`ADDL`, `MDBDS`, `MDTTE`,
`MDOCCDS`) belong to a device rather than a subject, so devices are
counted by `SPDEVID`, and the per-parameter counts and the checks are
per device. `USUBJID` is optional in these datasets; if it is present
the subjects are counted as well, and it is what `adsl` is matched by.

The dataset name is the one stored by
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
or, failing that, the name of the variable
[`summary()`](https://rdrr.io/r/base/summary.html) was called on.

## See also

Utilities for Dataset Checking:
[`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_duplicates_dataset.md),
[`get_many_to_one_dataset()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_many_to_one_dataset.md),
[`get_one_to_many_dataset()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_one_to_many_dataset.md)

## Examples

### Summarizing a BDS dataset after a derivation

    library(tibble)
    library(dplyr, warn.conflicts = FALSE)

    advs <- tribble(
      ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
      "01-701-1015", "DIABP",     51, "BASELINE",
      "01-701-1015", "SYSBP",    121, "BASELINE",
      "01-701-1028", "DIABP",     79, "BASELINE",
      "01-701-1028", "SYSBP",    130, "BASELINE"
    ) %>%
      mutate(STUDYID = "PILOT01")

    advs <- derive_param_computed(
      advs,
      by_vars = exprs(STUDYID, USUBJID, AVISIT),
      parameters = c("SYSBP", "DIABP"),
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP"
      )
    )

    summary(advs)
    #> ── advs: BDS summary ───────────────────────────────────────────────────────────
    #> Subjects (USUBJID): 2
    #> Observations: 6 | Variables: 5
    #> i Structure: not checked, no keys set (see `admiral_df()`)
    #> Parameters (3):
    #>   PARAMCD  records  subjects  visits  missing    min  median  max
    #>   DIABP          2         2       1        0     51      65   79
    #>   MAP            2         2       1        0  74.33   85.17   96
    #>   SYSBP          2         2       1        0    121   125.5  130
    #> Analysis visits (AVISIT): BASELINE

### Checking against declared keys

The record structure is only checked once the keys have been stored with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).

    advs <- admiral_df(
      advs,
      dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS"
    )

    summary(advs)
    #> ── ADVS: BDS summary ───────────────────────────────────────────────────────────
    #> Subjects (USUBJID): 2
    #> Observations: 6 | Variables: 5
    #> ✔ Structure: one record per USUBJID, PARAMCD, AVISIT
    #> Parameters (3):
    #>   PARAMCD  records  subjects  visits  missing    min  median  max
    #>   DIABP          2         2       1        0     51      65   79
    #>   MAP            2         2       1        0  74.33   85.17   96
    #>   SYSBP          2         2       1        0    121   125.5  130
    #> Analysis visits (AVISIT): BASELINE

### Comparing with ADSL

Supplying the ADSL shows the subjects without records, and checks that
every subject of the dataset is in the ADSL.

    adsl <- tribble(
      ~STUDYID,  ~USUBJID,      ~SAFFL, ~TRT01A,
      "PILOT01", "01-701-1015", "Y",    "Placebo",
      "PILOT01", "01-701-1028", "Y",    "Xanomeline High Dose",
      "PILOT01", "01-701-1034", "N",    NA
    )

    summary(advs, adsl = adsl)
    #> ── ADVS: BDS summary ───────────────────────────────────────────────────────────
    #> Subjects (USUBJID): 2
    #> Observations: 6 | Variables: 5
    #> ✔ Structure: one record per USUBJID, PARAMCD, AVISIT
    #> Parameters (3):
    #>   PARAMCD  records  subjects  visits  missing    min  median  max
    #>   DIABP          2         2       1        0     51      65   79
    #>   MAP            2         2       1        0  74.33   85.17   96
    #>   SYSBP          2         2       1        0    121   125.5  130
    #> Analysis visits (AVISIT): BASELINE
    #> Subjects with records out of ADSL: 2 of 3 (66.7%) | safety population 2 of 2
    #> (100.0%)
    #> Subjects with records per arm (TRT01A): Placebo: 1 of 1 | Xanomeline High Dose:
    #> 1 of 1 | NA: 0 of 1
    #> ✔ Checks passed: all subjects with records are in ADSL
