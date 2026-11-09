# admiral_df Test 9: empty keys are stored, keys not in the dataset are not

    Code
      typo <- admiral_df(input, dataset_keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
    Condition
      Warning:
      Key variable `PARMACD` is not in `dataset`, so no keys were set.

---

    Code
      rekeyed <- admiral_df(keyed, dataset_keys = exprs(USUBJID, PARMACD, AVISIT))
    Condition
      Warning:
      Key variables `PARMACD` and `AVISIT` are not in `dataset`, so no keys were set.
      i The keys set previously are kept: `USUBJID` and `PARAMCD`.

# admiral_df Test 11: an error is issued for invalid arguments

    Code
      admiral_df(input, dataset_class = "ADAE")
    Condition
      Error in `admiral_df()`:
      ! Argument `dataset_class` must be equal to one of "ADSL", "BDS", "OCCDS", "TTE", "ADDL", "MDBDS", "MDTTE", or "MDOCCDS".

# print.admiral_df Test 27: the name, class and keys are printed in a header

    Code
      print(admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS", dataset_class = "BDS"))
    Output
      # An admiral dataset (name: ADVS, class: BDS, keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT)))
    Output
      # An admiral dataset (class: BDS (inferred), keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, dataset_class = "BDS"))
    Output
      # An admiral dataset (class: BDS)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, dataset_keys = exprs(USUBJID), dataset_class = "TTE"))
    Output
      # An admiral dataset (class: TTE, keys: USUBJID)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(lookup, dataset_keys = exprs(PARAMCD), dataset_name = "PARAMS"))
    Output
      # An admiral dataset (name: PARAMS, keys: PARAMCD)
      # A tibble: 1 x 2
        PARAMCD PARAM                          
        <chr>   <chr>                          
      1 DIABP   Diastolic Blood Pressure (mmHg)

---

    Code
      print(admiral_df(input, dataset_keys = exprs(), dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, class: BDS (inferred), keys: none)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(as.data.frame(input), dataset_keys = exprs(USUBJID),
      dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, class: BDS (inferred), keys: USUBJID)
        STUDYID USUBJID PARAMCD   AVISIT AVAL
      1       P       1   DIABP BASELINE   51
      2       P       1   SYSBP BASELINE  121

# summary.admiral_df Test 30: the structure is not checked if no keys were set

    Code
      print(actual)
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          3         2       2        0  118     121  130
      Analysis visits (AVISIT): BASELINE, WEEK 2

# summary.admiral_df Test 31: duplicate, missing and empty declared keys

    Code
      print(duplicated)
    Message
      -- BDS summary -----------------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
      x Structure: expected one record per USUBJID, PARAMCD, found 1 extra record
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          3         2       2        0  118     121  130
      Analysis visits (AVISIT): BASELINE, WEEK 2

---

    Code
      print(summary(renamed))
    Message
      -- renamed: BDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
      ! Structure: not checked, key AVISIT is not in the dataset
      Parameters (1):
        PARAMCD  records  subjects  missing  min  median  max
        SYSBP          3         2        0  118     121  130

---

    Code
      print(summary(renamed_two))
    Message
      -- renamed_two: Dataset summary ------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
      ! Structure: not checked, keys PARAMCD, AVISIT are not in the dataset

---

    Code
      print(no_keys)
    Message
      -- BDS summary -----------------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
      i Structure: no key variables defined
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          3         2       2        0  118     121  130
      Analysis visits (AVISIT): BASELINE, WEEK 2

# summary.admiral_df Test 32: the ADSL summary

    Code
      print(actual)
    Message
      -- ADSL summary ----------------------------------------------------------------
      Subjects (USUBJID): 3
      Observations: 3 | Variables: 5
      v Structure: one record per USUBJID
      Subjects per treatment arm (TRT01P): Drug A: 2 | Placebo: 1
      Populations: SAFFL 2 | ITTFL 3

# summary.admiral_df Test 33: the BDS summary

    Code
      print(actual)
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 4 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          1         1       1        0   79      79   79
        SYSBP          3         2       2        0  118     121  130
      Analysis visits (AVISIT): BASELINE, WEEK 2

---

    Code
      print(many_summary)
    Message
      -- many_params: BDS summary ----------------------------------------------------
      Subjects (USUBJID): 1
      Observations: 12 | Variables: 4
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (12):
        PARAMCD  records  subjects  missing  min  median  max
        P01            1         1        0    1       1    1
        P02            1         1        0    1       1    1
        P03            1         1        0    1       1    1
        P04            1         1        0    1       1    1
        P05            1         1        0    1       1    1
        P06            1         1        0    1       1    1
        P07            1         1        0    1       1    1
        P08            1         1        0    1       1    1
        P09            1         1        0    1       1    1
        P10            1         1        0    1       1    1
      ... and 2 more parameters

# summary.admiral_df Test 34: the OCCDS summary

    Code
      print(actual)
    Message
      -- input: OCCDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Distinct terms: AEBODSYS 2 | AEDECOD 2
      Treatment emergent (TRTEMFL): 2 of 3 records

# summary.admiral_df Test 35: the TTE summary

    Code
      print(actual)
    Message
      -- input: TTE summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Events and censored per parameter (CNSR):
        PARAMCD  events  censored
        OS            1         1
        PFS           0         1

# summary.admiral_df Test 36: datasets with little to report do not error

    Code
      print(summary(lookup))
    Message
      -- lookup: Dataset summary -----------------------------------------------------
      Observations: 2 | Variables: 1
      i Structure: not checked, no keys set (see `admiral_df()`)

---

    Code
      print(summary(empty))
    Message
      -- empty: Dataset summary ------------------------------------------------------
      Subjects (USUBJID): 0
      Observations: 0 | Variables: 2
      i Structure: not checked, no keys set (see `admiral_df()`)

# print_summary_table Test 38: tables are aligned, counts are not rounded, NA is shown as `.`

    Code
      print_summary_table(tbl)
    Message
        PARAMCD  records   mean  subjects
        DIABP      11509  80.12       254
        SYSBP          7      .         3
        .              1      1         .

---

    Code
      print_summary_table(tbl, max_rows = 1, what = "parameter")
    Message
        PARAMCD  records   mean  subjects
        DIABP      11509  80.12       254
      ... and 2 more parameters

# summary.admiral_df Test 39: the ADSL subject flow, outcome and demographics

    Code
      print(summary(input))
    Message
      -- ADSL summary ----------------------------------------------------------------
      Subjects (USUBJID): 4
      Observations: 4 | Variables: 10
      v Structure: one record per USUBJID
      Subject flow: randomized 3 | treated 2 | completed 1
      Treatment duration (TRTDURD, days): median 60 (range 20-100)
      End of study (EOSSTT): COMPLETED: 1 | DISCONTINUED: 2 | NA: 1
      Deaths (DTHFL): 1 (ADVERSE EVENT: 1)
      Demographics: AGE median 59 (range 40-71) | SEX F: 2 | M: 2

# summary.admiral_df Test 40: the BDS per-parameter statistics and derived records

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 5 | Variables: 7
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median    max
        DIABP          1         1       0        1    .       .      .
        SYSBP          4         2       2        1  120     125  130.5
      Analysis visits (AVISIT): BASELINE, WEEK 2
      Derived records (DTYPE): LOCF: 1

# summary.admiral_df Test 42: more than one baseline per subject and parameter is reported

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 5 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          2         1       2        0   80      81   82
        SYSBP          3         2       3        0  120     121  130
      Analysis visits (AVISIT): BASELINE, SCREENING, WEEK 2
      x 1 combination of USUBJID, PARAMCD has more than one baseline record (`ABLFL`)

---

    Code
      print(summary(input[-1, ]))
    Message
      -- BDS summary -----------------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 4 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          2         1       2        0   80      81   82
        SYSBP          2         2       2        0  121   125.5  130
      Analysis visits (AVISIT): BASELINE, WEEK 2
      v Checks passed: at most one baseline record per USUBJID, PARAMCD

# print_summary_checks Test 43: failed checks are listed, passed ones combined, absent ones left out

    Code
      print_summary_checks(d, checks)
    Message
      x 2 a failures
      v Checks passed: no b; no c

# summary.admiral_df Test 44: the OCCDS severity and seriousness

    Code
      print(summary(input))
    Message
      -- input: OCCDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 7
      i Structure: not checked, no keys set (see `admiral_df()`)
      Distinct terms: AEDECOD 2
      Records per severity (ASEV): MILD: 1 | SEVERE: 1 | NA: 1
      Serious (AESER): 1 of 3 records

# check_occurrence_flags Test 45: each flag is checked within the group its name defines

    Code
      print(summary(admiral_df(dplyr::mutate(dup, STUDYID = "P"))))
    Message
      -- OCCDS summary ---------------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 6 | Variables: 9
      i Structure: not checked, no keys set (see `admiral_df()`)
      Distinct terms: AEBODSYS 2 | AEDECOD 3
      x 4 occurrence flags have more than one "Y" per group:
        AOCCFL per USUBJID: 1 group
        AOCCSFL per USUBJID, AEBODSYS: 1 group
        AOCCPFL per USUBJID, AEDECOD: 2 groups
        AOCCPIFL per USUBJID, AEDECOD: 2 groups
      i Not checked, as their level is not known: `AOCC01FL`

# summary.admiral_df Test 46: the dataset is compared with the supplied ADSL

    Code
      print(summary(adae, adsl = adsl))
    Message
      -- adae: OCCDS summary ---------------------------------------------------------
      Subjects (USUBJID): 3
      Observations: 5 | Variables: 4
      i Structure: not checked, no keys set (see `admiral_df()`)
      Distinct terms: AEDECOD 2
      Treatment emergent (TRTEMFL): 4 of 5 records
      Subjects with records out of ADSL: 2 of 4 (50.0%) | safety population 2 of 3
      (66.7%)
      Subjects with records per arm (TRT01A): Drug A: 1 of 2 | Placebo: 1 of 1 | NA:
      0 of 1
      Subjects with a treatment-emergent record (TRTEMFL): 1 of 3 subjects in the
      safety population (33.3%)
      x 1 subject with records is not in ADSL: "9"

---

    Code
      print(summary(advs, adsl = adsl))
    Message
      -- advs: BDS summary -----------------------------------------------------------
      Subjects (USUBJID): 1
      Observations: 1 | Variables: 4
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (1):
        PARAMCD  records  subjects  missing  min  median  max
        SYSBP          1         1        0  120     120  120
      Subjects with records out of ADSL: 1 of 4 (25.0%) | safety population 1 of 3
      (33.3%)
      Subjects with records per arm (TRT01A): Drug A: 0 of 2 | Placebo: 1 of 1 | NA:
      0 of 1
      v Checks passed: all subjects with records are in ADSL

---

    Code
      print(summary(many_orphans, adsl = adsl))
    Message
      -- many_orphans: OCCDS summary -------------------------------------------------
      Subjects (USUBJID): 10
      Observations: 10 | Variables: 3
      i Structure: not checked, no keys set (see `admiral_df()`)
      Distinct terms: AEDECOD 1
      Subjects with records out of ADSL: 0 of 4 (0.0%) | safety population 0 of 3
      (0.0%)
      Subjects with records per arm (TRT01A): Drug A: 0 of 2 | Placebo: 0 of 1 | NA:
      0 of 1
      x 10 subjects with records are not in ADSL: "11", "12", "13", ..., "19", and
        "20"

# summary.admiral_df Test 47: an ADSL which cannot be compared with is an error

    Code
      summary(adae, adsl = adsl[c(1, 1, 2), ])
    Condition
      Error in `summary()`:
      ! `adsl` must have one record per subject (`USUBJID`).
      i Found 1 extra subject record.

---

    Code
      actual <- summary(admiral_df(adsl), adsl = adsl)
    Condition
      Warning:
      `adsl` is ignored for a dataset of class "ADSL", as it is meant for comparing another dataset with its ADSL.

# summary.admiral_df Test 50: a device dataset is summarized per device

    Code
      print(actual)
    Message
      -- mdbds: MDBDS summary --------------------------------------------------------
      Devices (SPDEVID): 2
      Observations: 4 | Variables: 6
      i Structure: not checked, no keys set (see `admiral_df()`)
      Parameters (1):
        PARAMCD  records  devices  visits  missing  min  median  max
        BATTV          4        2       2        0  2.8    2.95  3.1
      Analysis visits (AVISIT): BASELINE, WEEK 4
      v Checks passed: at most one baseline record per SPDEVID, PARAMCD

---

    Code
      print(addl_actual)
    Message
      -- addl: ADDL summary ----------------------------------------------------------
      Devices (SPDEVID): 3
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 4
      i Structure: not checked, no keys set (see `admiral_df()`)
      Subjects with records out of ADSL: 2 of 3 (66.7%)
      v Checks passed: all subjects with records are in ADSL

