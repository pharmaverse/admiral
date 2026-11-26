# Summarizing ADaM Datasets

## Introduction

While an ADaM dataset is being derived it is easy to lose track of what
it contains: did
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/derive_param_computed.md)
add the parameter you expected, is the correct number of flags present
in the dataset, is there still only one baseline record per subject and
parameter?

[admiral](https://pharmaverse.github.io/admiral/) helps answer these
questions with the `admiral_df` class. A dataset of this class can carry
a little metadata about itself – its name, its ADaM class, and its key
variables – and calling
[`summary()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/summary.admiral_df.md)
on it prints a short report tailored to the kind of dataset it is,
including some sanity checks.

This article shows how to:

- Use [`summary()`](https://rdrr.io/r/base/summary.html) to inspect a
  dataset *in flight*, i.e. at any point of a derivation pipeline,
- Tag a dataset with
  [`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
  and read its metadata back,
- Keep the class and metadata when using
  [dplyr](https://dplyr.tidyverse.org) and
  [tidyr](https://tidyr.tidyverse.org), as some common operations keep
  them, some keep them only in part, and some drop them.

**Note**: The `admiral_df` class and the functions described here are
[experimental](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#experimental-functions),
so their behavior may still change.

### Required Packages

The examples in this article use the following packages and the example
ADSL and ADLB datasets which come with
[admiral](https://pharmaverse.github.io/admiral/). Both datasets already
have their metadata stored, so they can be summarized straight away.

``` r
library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(tibble)

adsl <- admiral_adsl
adlb <- admiral_adlb
```

## Summarizing a Dataset in Flight

The example datasets already have their name, ADaM class, and key
variables stored (see [Tagging a Dataset](#tagging-a-dataset)), so
[`summary()`](https://rdrr.io/r/base/summary.html) can report on them
straight away.

Calling [`summary()`](https://rdrr.io/r/base/summary.html) on an
`admiral_df` prints a report on the dataset. It always shows the name
and ADaM class, the number of subjects, records and variables, and
whether there is one record per key. Depending on the ADaM class it then
adds:

- **ADSL**: subjects per treatment arm, population flags, subject flow
  (randomized, treated, completed), treatment duration, end of study
  status, deaths and their causes, and age and sex.
- **BDS**: the records, subjects and visits of each parameter, with the
  distribution of `AVAL` for numeric parameters and the distinct `AVALC`
  values for character ones, the analysis visits, the derived records
  per `DTYPE`, and checks that every parameter has `AVAL` or `AVALC`
  populated and that there is at most one baseline record (`ABLFL`) per
  subject, parameter, and `BASETYPE`.
- **OCCDS**: the number of distinct terms, treatment-emergent and
  serious records, records per severity and per toxicity grade, and a
  check that each occurrence flag (`AOCCFL`, `AOCCPFL`, …) has at most
  one `"Y"` per group.
- **TTE**: events and censored records per parameter.

``` r
summary(adsl)
#> ── ADSL summary ────────────────────────────────────────────────────────────────
#> Subjects (USUBJID): 306
#> Observations: 306 | Variables: 54
#> ✔ Structure: one record per STUDYID, USUBJID
#> Subjects per treatment arm (TRT01P): Placebo: 86 | Screen Failure: 52 |
#> Xanomeline High Dose: 84 | Xanomeline Low Dose: 84
#> Populations: SAFFL 254
#> Subject flow: randomized 254 | treated 254 | completed 110
#> Treatment duration (TRTDURD, days): median 132 (range 1-212)
#> End of study (EOSSTT): COMPLETED: 110 | DISCONTINUED: 144 | NA: 52
#> Deaths (DTHFL): 3 (COMPLETED SUICIDE: 1 | MYOCARDIAL INFARCTION: 1 | SUDDEN
#> DEATH: 1)
#> Demographics: AGE median 77 (range 50-89) | SEX F: 179 | M: 127
```

``` r
summary(adlb)
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Subjects (USUBJID): 11
#> Observations: 330 | Variables: 115
#> ✔ Structure: one record per STUDYID, USUBJID, PARAMCD, ADT, AVISITN, VISITNUM
#> 3 numeric parameters (AVAL):
#>   PARAMCD  records  subjects  visits  missing   min  median    max
#>   ALT          110        11      18        0    10      23    129
#>   AST          110        11      18        0    16      24    168
#>   BILI         110        11      18        0  5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✔ Checks passed: every parameter has AVAL or AVALC populated; at most one
#>   baseline record per USUBJID, PARAMCD, BASETYPE
```

Each item is only reported when the variables it needs are in the
dataset. So [`summary()`](https://rdrr.io/r/base/summary.html) can be
called at any point of a derivation: a dataset which is only half
derived gives a shorter report, not an error.

### Checking the Record Structure

The structure line compares the records with the key variables stored
with the dataset:

- `✔` – there is one record per combination of the keys,
- `✖` – some combinations have more than one record, and how many extra
  records there are,
- `ℹ` – no keys have been set, so the structure was not checked,
- `!` – some keys are no longer variables of the dataset (e.g. after a
  [`select()`](https://dplyr.tidyverse.org/reference/select.html) or
  [`rename()`](https://dplyr.tidyverse.org/reference/rename.html)), so
  the structure was not checked.

For example, imagine a merge accidentally duplicated some records.
Without looking through the data,
[`summary()`](https://rdrr.io/r/base/summary.html) reports it straight
away:

``` r
adlb_dup <- bind_rows(adlb, slice(adlb, 1:2))

summary(adlb_dup)
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Subjects (USUBJID): 11
#> Observations: 332 | Variables: 115
#> ✖ Structure: expected one record per STUDYID, USUBJID, PARAMCD, ADT, AVISITN,
#>   VISITNUM, found 2 extra records
#> 3 numeric parameters (AVAL):
#>   PARAMCD  records  subjects  visits  missing   min  median    max
#>   ALT          112        11      18        0    10      23    129
#>   AST          110        11      18        0    16      24    168
#>   BILI         110        11      18        0  5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✖ 1 combination of USUBJID, PARAMCD, BASETYPE has more than one baseline record
#>   (`ABLFL`)
#> ✔ Checks passed: every parameter has AVAL or AVALC populated
```

### Comparing with ADSL

A dataset cannot tell by itself which subjects are missing from it. If
ADSL is passed to the `adsl` argument,
[`summary()`](https://rdrr.io/r/base/summary.html) also reports how many
of the ADSL subjects have records – overall, in the safety population,
and per treatment arm – and checks that every subject in the dataset is
in ADSL. For an OCCDS dataset it also reports the subjects with a
treatment-emergent record.

``` r
summary(adlb, adsl = adsl)
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Subjects (USUBJID): 11
#> Observations: 330 | Variables: 115
#> ✔ Structure: one record per STUDYID, USUBJID, PARAMCD, ADT, AVISITN, VISITNUM
#> 3 numeric parameters (AVAL):
#>   PARAMCD  records  subjects  visits  missing   min  median    max
#>   ALT          110        11      18        0    10      23    129
#>   AST          110        11      18        0    16      24    168
#>   BILI         110        11      18        0  5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✔ Checks passed: every parameter has AVAL or AVALC populated; at most one
#>   baseline record per USUBJID, PARAMCD, BASETYPE
#> Subjects with records out of ADSL: 11 of 306 (3.6%) | safety population 11 of
#> 254 (4.3%)
#> Subjects with records per arm (TRT01A): Placebo: 5 of 86 | Screen Failure: 0 of
#> 52 | Xanomeline High Dose: 3 of 72 | Xanomeline Low Dose: 3 of 96
#> ✔ Checks passed: all subjects with records are in ADSL
```

### A Worked Example

Summaries are most useful straight after a derivation, to confirm it did
what you expected. Here we add a parameter for the ratio of ALT to AST.
The keys are set to one record per subject, parameter, date, and `DTYPE`
first, so that the derived records (`DTYPE` is populated) and the
observed ones are told apart (see [Tagging a
Dataset](#tagging-a-dataset)):

``` r
adlb <- admiral_df(adlb, dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, ADT, DTYPE))

adlb_ratio <- adlb %>%
  derive_param_computed(
    by_vars = exprs(STUDYID, USUBJID, AVISIT, AVISITN, ADT),
    parameters = c("ALT", "AST"),
    set_values_to = exprs(
      AVAL = AVAL.ALT / AVAL.AST,
      PARAMCD = "ALTAST",
      PARAM = "ALT/AST Ratio"
    )
  )

summary(adlb_ratio)
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Subjects (USUBJID): 11
#> Observations: 436 | Variables: 115
#> ✖ Structure: expected one record per STUDYID, USUBJID, PARAMCD, ADT, DTYPE,
#>   found 26 extra records
#> 4 numeric parameters (AVAL):
#>   PARAMCD  records  subjects  visits  missing     min  median    max
#>   ALT          110        11      18        0      10      23    129
#>   ALTAST       106        11      18        0  0.4444  0.8297      2
#>   AST          110        11      18        0      16      24    168
#>   BILI         110        11      18        0    5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✔ Checks passed: every parameter has AVAL or AVALC populated; at most one
#>   baseline record per USUBJID, PARAMCD, BASETYPE
```

The derivation ran without an error or a warning, but the summary shows
that something went wrong: the structure check finds 26 extra records.
The ratio was also computed for the derived `"POST-BASELINE"` visits
(`DTYPE` is `"LOV"`, `"MAXIMUM"`, or `"MINIMUM"`). As `DTYPE` is not set
for the new records, these duplicate the ratio of the observed visit on
the same date.

The fix is to compute the ratio from the observed records only. The
`filter` argument of
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/derive_param_computed.md)
does this without dropping the derived records from the dataset, as
filtering the dataset beforehand would:

``` r
adlb_ratio <- adlb %>%
  derive_param_computed(
    by_vars = exprs(STUDYID, USUBJID, AVISIT, AVISITN, ADT),
    filter = is.na(DTYPE),
    parameters = c("ALT", "AST"),
    set_values_to = exprs(
      AVAL = AVAL.ALT / AVAL.AST,
      PARAMCD = "ALTAST",
      PARAM = "ALT/AST Ratio"
    )
  )

summary(adlb_ratio)
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Subjects (USUBJID): 11
#> Observations: 410 | Variables: 115
#> ✔ Structure: one record per STUDYID, USUBJID, PARAMCD, ADT, DTYPE
#> 4 numeric parameters (AVAL):
#>   PARAMCD  records  subjects  visits  missing     min  median    max
#>   ALT          110        11      18        0      10      23    129
#>   ALTAST        80        11      15        0  0.4444  0.8367      2
#>   AST          110        11      18        0      16      24    168
#>   BILI         110        11      18        0    5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✔ Checks passed: every parameter has AVAL or AVALC populated; at most one
#>   baseline record per USUBJID, PARAMCD, BASETYPE
```

Now the structure check passes, and `ALTAST` has a record for each
observed visit with both values.

### Using the Summary Programmatically

[`summary()`](https://rdrr.io/r/base/summary.html) returns a list, which
is printed. Its elements can be used in your own checks, e.g. to warn,
or stop a program, if the record structure is wrong:

``` r
smry <- summary(adlb_dup)

smry$structure$n_duplicates
#> [1] 2

if (isTRUE(smry$structure$n_duplicates > 0)) {
  cli::cli_warn("Duplicates detected!")
}
#> Warning: Duplicates detected!
```

`n_duplicates` is `NA` if the structure was not checked, e.g. because no
keys were set, hence the
[`isTRUE()`](https://rdrr.io/r/base/Logic.html).

## What is an `admiral_df`?

The summaries above rely on the `admiral_df` class and the metadata
stored with it. An `admiral_df` is an ordinary data frame (usually a
tibble) with an extra class. Everything you can do with a tibble you can
still do with an `admiral_df`; it only changes how the dataset is
printed and summarized.

All [admiral](https://pharmaverse.github.io/admiral/) functions which
create or modify a dataset, such as the `derive_*()` functions, return
an `admiral_df`. So as soon as a dataset has passed through an
[admiral](https://pharmaverse.github.io/admiral/) function,
[`summary()`](https://rdrr.io/r/base/summary.html) gives the
[admiral](https://pharmaverse.github.io/admiral/) report rather than the
base R one.

``` r
adsl <- derive_vars_dt(
  adsl,
  new_vars_prefix = "RFST",
  dtc = RFSTDTC
)

class(adsl)
#> [1] "admiral_df" "tbl_df"     "tbl"        "data.frame"
```

In addition, a dataset can store three pieces of metadata:

| Metadata | Example | Read back with |
|----|----|----|
| Name | `"ADLB"` | [`get_dataset_name()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_name.md) |
| ADaM class | `"BDS"` | [`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_class.md) |
| Key variables | `exprs(USUBJID, PARAMCD, AVISIT)` | [`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_keys.md) |

The key variables are the variables which together should uniquely
identify each record of the dataset. They are what
[`summary()`](https://rdrr.io/r/base/summary.html) checks the record
structure against.

## Tagging a Dataset

The metadata is set with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).
Any data frame can be passed to it, whether or not it has the
`admiral_df` class yet.

``` r
adlb <- admiral_df(
  adlb,
  dataset_name = "ADLB",
  dataset_class = "BDS",
  dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, AVISIT)
)
```

When the metadata has been set, it is shown in a header line when the
dataset is printed:

``` r
adlb %>%
  select(USUBJID, PARAMCD, AVISIT, AVAL) %>%
  head(3)
#> # An admiral dataset (name: ADLB, class: BDS, keys: STUDYID, USUBJID, PARAMCD, AVISIT)
#> # A tibble: 3 × 4
#>   USUBJID     PARAMCD AVISIT                 AVAL
#>   <chr>       <chr>   <chr>                 <dbl>
#> 1 01-701-1015 ALT     Baseline                 27
#> 2 01-701-1015 ALT     Week 2                   41
#> 3 01-701-1015 ALT     POST-BASELINE MAXIMUM    41
```

and it can be read back with the `get_*()` functions listed in the table
above.
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/get_dataset_keys.md)
returns a list of symbols, so the keys can be passed straight to a
`by_vars` argument of another function:

``` r
get_dataset_keys(adlb)
#> [[1]]
#> STUDYID
#> 
#> [[2]]
#> USUBJID
#> 
#> [[3]]
#> PARAMCD
#> 
#> [[4]]
#> AVISIT
```

### Updating the Metadata

Calling
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
again on a dataset which already has metadata updates only the items
which are supplied. Anything left as `NULL` is kept:

``` r
adlb <- admiral_df(adlb, dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, ADT, DTYPE))

get_dataset_name(adlb)
#> [1] "ADLB"
get_dataset_keys(adlb)
#> [[1]]
#> STUDYID
#> 
#> [[2]]
#> USUBJID
#> 
#> [[3]]
#> PARAMCD
#> 
#> [[4]]
#> ADT
#> 
#> [[5]]
#> DTYPE
```

If one of the keys is not a variable of the dataset, perhaps because of
a typo, a warning is issued and the keys are left as they were:

``` r
adlb_typo <- admiral_df(adlb, dataset_keys = exprs(USUBJID, PARAMDC))
#> Warning: Key variable `PARAMDC` is not in `dataset`, so no keys were set.
#> ℹ The keys set previously are kept: `STUDYID`, `USUBJID`, `PARAMCD`, `ADT`, and
#>   `DTYPE`.
```

To record that a dataset has no key variables at all, pass an empty
list,
[`exprs()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/reexport-exprs.md).
This is different from never setting the keys.

### Where do the Keys Come From?

[admiral](https://pharmaverse.github.io/admiral/) never guesses the keys
of a dataset. They can come from your dataset specification, from the
`by_vars` of the derivation which created the dataset, or simply be
written by hand. A specification usually provides them as a character
vector, which
[`chr2vars()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/chr2vars.md)
converts. For example, with a
[metacore](https://atorus-research.github.io/metacore/) specification:

``` r
adlb <- admiral_df(
  adlb,
  dataset_keys = chr2vars(metacore::get_keys(spec, "ADLB")$variable),
  dataset_name = "ADLB"
)
```

### The ADaM Class

The ADaM class decides what
[`summary()`](https://rdrr.io/r/base/summary.html) reports. The
supported classes are `"ADSL"`, `"BDS"`, `"OCCDS"`, and `"TTE"`, and for
medical device data `"ADDL"`, `"MDBDS"`, `"MDOCCDS"`, and `"MDTTE"`.

Unlike the keys, the class does not have to be set: if it has not been,
[admiral](https://pharmaverse.github.io/admiral/) does its best to infer
it from the variables and records of the dataset. For example, a dataset
with `PARAMCD` and `AVAL` is inferred as `"BDS"`, one with occurrence
variables such as `AEDECOD` or `TRTEMFL` as `"OCCDS"`, and one with
neither and one record per subject as `"ADSL"`:

``` r
adae <- tribble(
  ~STUDYID,  ~USUBJID,      ~AEDECOD,
  "PILOT01", "01-701-1015", "HEADACHE",
  "PILOT01", "01-701-1023", "NAUSEA"
)

get_dataset_class(adae)
#> [1] "OCCDS"
```

When printed, an inferred class is marked as such. If the inferred class
is wrong for your dataset – for example an ADAE which has been
restricted to device related events may be taken for device data – set
it with
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md).
A class that has been set always takes precedence.

## What Happens to the Class and Metadata?

The class and metadata are stored as attributes of the data frame.
Whether they survive an operation depends on whether it modifies the
existing data frame or builds a new one. The table below summarizes what
common operations do to an `admiral_df` with metadata set.

| Operation | Class | Metadata |
|:---|:---|:---|
| [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html), [`filter()`](https://dplyr.tidyverse.org/reference/filter.html), [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html), [`slice()`](https://dplyr.tidyverse.org/reference/slice.html), [`relocate()`](https://dplyr.tidyverse.org/reference/relocate.html), [`head()`](https://rdrr.io/r/utils/head.html) | kept | kept |
| `mutate(.by = )`, `filter(.by = )` | kept | kept |
| [`left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html) etc. and [`semi_join()`](https://dplyr.tidyverse.org/reference/filter-joins.html)/[`anti_join()`](https://dplyr.tidyverse.org/reference/filter-joins.html), dataset as `x` | kept | kept |
| [`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html), dataset as *first* argument; [`rbind()`](https://rdrr.io/r/base/cbind.html) | kept | kept |
| [`select()`](https://dplyr.tidyverse.org/reference/select.html), [`rename()`](https://dplyr.tidyverse.org/reference/rename.html), [`transmute()`](https://dplyr.tidyverse.org/reference/transmute.html), [`distinct()`](https://dplyr.tidyverse.org/reference/distinct.html), [`count()`](https://dplyr.tidyverse.org/reference/count.html) | kept | kept, may be out of date |
| [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html), [`rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html) | dropped | kept |
| [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html), [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) | dropped | kept |
| [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) … [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html), [`rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html) … [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html) | dropped | dropped |
| [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html), [`reframe()`](https://dplyr.tidyverse.org/reference/reframe.html) | dropped | dropped |
| [`left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html) etc., dataset as `y` | dropped | dropped |
| [`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html), dataset *not* the first argument | dropped | dropped |
| [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html), [`pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html), [`complete()`](https://tidyr.tidyverse.org/reference/complete.html), [`nest()`](https://tidyr.tidyverse.org/reference/nest.html)/[`unnest()`](https://tidyr.tidyverse.org/reference/unnest.html) | dropped | dropped |
| [admiral](https://pharmaverse.github.io/admiral/) functions, e.g. [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/derive_var_extreme_flag.md) | kept | not guaranteed |

### Kept

Operations which change the values or the order of the records, but keep
the data frame they were given, keep everything. This covers most of the
code in a typical ADaM script:

``` r
adlb %>%
  filter(PARAMCD == "ALT") %>%
  mutate(AVALCAT1 = if_else(AVAL > 40, "> 40", "<= 40")) %>%
  get_dataset_name()
#> [1] "ADLB"
```

Joins and
[`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
keep the attributes of the dataset which comes *first*. If the dataset
with the metadata is the `x` of a join, or the first dataset of
[`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html),
everything is kept; otherwise the metadata of the other dataset is used,
which usually means none.

``` r
# ADLB is `x`: kept
left_join(adlb, select(adsl, STUDYID, USUBJID, RFSTDT), by = c("STUDYID", "USUBJID")) %>%
  get_dataset_name()
#> [1] "ADLB"

# ADLB is `y`: dropped
left_join(select(adsl, STUDYID, USUBJID), adlb, by = c("STUDYID", "USUBJID")) %>%
  get_dataset_name()
#> [1] "ADSL"
```

### Partially Kept

#### Metadata kept, but possibly out of date

[`select()`](https://dplyr.tidyverse.org/reference/select.html),
[`rename()`](https://dplyr.tidyverse.org/reference/rename.html),
[`transmute()`](https://dplyr.tidyverse.org/reference/transmute.html),
[`distinct()`](https://dplyr.tidyverse.org/reference/distinct.html), and
[`count()`](https://dplyr.tidyverse.org/reference/count.html) keep the
class and metadata, but they change the variables, so the stored keys
may no longer match the dataset. The keys are stored as plain names, so
they are not renamed along with the variables.
[`summary()`](https://rdrr.io/r/base/summary.html) tells you when this
has happened:

``` r
adlb %>%
  rename(SUBJECT = USUBJID) %>%
  summary()
#> ── ADLB: BDS summary ───────────────────────────────────────────────────────────
#> Observations: 330 | Variables: 115
#> ! Structure: not checked, key USUBJID is not in the dataset
#> 3 numeric parameters (AVAL):
#>   PARAMCD  records  visits  missing   min  median    max
#>   ALT          110      18        0    10      23    129
#>   AST          110      18        0    16      24    168
#>   BILI         110      18        0  5.13    8.55  124.8
#> Analysis visits (AVISIT): Baseline, Week 2, Unscheduled 4.1, Unscheduled 4.2,
#> Week 4, Unscheduled 5.1, Unscheduled 6.1, Week 6, Week 8, Unscheduled 8.2, Week
#> 12, Week 16, Week 20, Week 24, Week 26, POST-BASELINE MINIMUM, POST-BASELINE
#> MAXIMUM, POST-BASELINE LAST
#> Derived records (DTYPE): LOV: 30 | MAXIMUM: 30 | MINIMUM: 30
#> ✔ Checks passed: every parameter has AVAL or AVALC populated
```

Similarly, [`count()`](https://dplyr.tidyverse.org/reference/count.html)
and [`distinct()`](https://dplyr.tidyverse.org/reference/distinct.html)
return a dataset which keeps the name and class of the original, even
though it is no longer an ADLB. Bear this in mind when printing or
summarizing their results.

#### Metadata kept, class dropped

[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) and
[`rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html) return
a grouped data frame, and
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
return a plain one. They keep the metadata, but not the `admiral_df`
class, so [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) fall back to their
usual behavior. The metadata is still there, so calling
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
without any other arguments restores the class together with it:

``` r
adlb_tbl <- as_tibble(adlb)

class(adlb_tbl)
#> [1] "tbl_df"     "tbl"        "data.frame"

adlb_tbl %>%
  admiral_df() %>%
  get_dataset_name()
#> [1] "ADLB"
```

#### `{admiral}` functions

All [admiral](https://pharmaverse.github.io/admiral/) functions return
an `admiral_df`, so the class is always kept. The metadata, however, is
not guaranteed to be: functions which build their result from a grouped
or summarized dataset internally lose it, while others keep it. For
example,
[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/derive_var_extreme_flag.md)
drops it:

``` r
adlb %>%
  filter(is.na(DTYPE)) %>%
  derive_var_extreme_flag(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    order = exprs(ADT),
    new_var = LASTFL,
    mode = "last"
  ) %>%
  get_dataset_name()
#> NULL
```

### Dropped

[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html),
the [tidyr](https://tidyr.tidyverse.org) reshaping functions, and the
grouped pipelines which end in
[`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html) build
a new data frame, which has neither the class nor the metadata. Note
that [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
followed by
[`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html) drops
everything, even though
[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) on
its own keeps the metadata:

``` r
adlb %>%
  group_by(USUBJID, PARAMCD) %>%
  mutate(N_VISITS = n()) %>%
  ungroup() %>%
  class()
#> [1] "tbl_df"     "tbl"        "data.frame"
```

Where possible, use the `.by` argument instead, which keeps everything:

``` r
adlb %>%
  mutate(N_VISITS = n(), .by = c(USUBJID, PARAMCD)) %>%
  get_dataset_name()
#> [1] "ADLB"
```

Otherwise, call
[`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
again with the metadata after the operation.

## Recommendations

- **Tag the dataset late.** As the metadata may be lost during a
  derivation, call
  [`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
  once the dataset is (nearly) finished, or just before you want to
  summarize it, rather than at the start of the script.
- **Re-tag after operations which drop it.** If you need the summary in
  the middle of a pipeline, call
  [`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
  again with the metadata. Keeping the keys in a variable,
  e.g. `adlb_keys <- exprs(STUDYID, USUBJID, PARAMCD, AVISIT)`, makes
  this easy.
- **Remove grouping first.**
  [`admiral_df()`](https:/pharmaverse.github.io/admiral/3209-5-summary-override/reference/admiral_df.md)
  does not accept grouped datasets, so call
  [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html)
  before it.
- **Summarize after each important derivation.** Comparing the summaries
  before and after a derivation is a quick way to see exactly what it
  added or changed.
- **Pass ADSL** to [`summary()`](https://rdrr.io/r/base/summary.html) to
  see the subjects which have no records, or should not be in the
  dataset.
