# summary.admiral_df Test 8: formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3
      v Structure (inferred): one record per USUBJID, PARAMCD
      Analysis visits (AVISIT): BASELINE
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          2         2       1        0   51      65   79
        SYSBP          1         1       1        0  121     121  121

# print.admiral_df_check Test 13: formatted output is stable

    Code
      print(check_admiral_df(input))
    Message
      -- admiral_df tag check --------------------------------------------------------
        is an admiral_df:   TRUE
        key variables:      USUBJID, PARAMCD
        keys present:       all present
        structure checked:  yes -- one record per USUBJID, PARAMCD (declared)

---

    Code
      print(check_admiral_df(dplyr::rename(input, SUBJID = USUBJID)))
    Message
      -- admiral_df tag check --------------------------------------------------------
        is an admiral_df:   TRUE
        key variables:      USUBJID, PARAMCD
        keys present:       stale: USUBJID missing; checked on PARAMCD only
        structure checked:  yes -- one record per PARAMCD (declared)

---

    Code
      print(check_admiral_df(tibble::as_tibble(input)))
    Message
      -- admiral_df tag check --------------------------------------------------------
        is an admiral_df:   FALSE
        key variables:      USUBJID, PARAMCD
        keys present:       all present
        structure checked:  no
      ! The `admiral_keys` attribute is orphaned: `summary()` dispatches to
        `summary.data.frame()`, so the diagnostic is lost. Re-tag with
        `as_admiral_df()`, passing `keys` explicitly to replace the existing
        attribute.

# print.summary_admiral_df Test 15: the subject-level check lines name what ran

    Code
      print(summary(clean))
    Message
      -- clean: ADSL summary ---------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 2
      v Structure (inferred): one record per USUBJID
      Treatment (ARM): Drug A: 1 | Placebo: 1
      Populations: SAFFL 2
      Subject flow: treated 2
      v Checks passed: ARM matches ACTARM; no missing TRTSDT in the safety
        population; population flags only contain Y/N/NA

---

    Code
      print(summary(broken))
    Message
      -- broken: ADSL summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 2
      v Structure (inferred): one record per USUBJID
      Treatment (ARM): Drug A: 1 | Placebo: 1
      Populations: SAFFL 2
      Subject flow: treated 1
      x 1 subject has different `ARM` and `ACTARM`
      x 1 subject in the safety population has no `TRTSDT`
      v Checks passed: population flags only contain Y/N/NA

---

    Code
      print(summary(mixed))
    Message
      -- mixed: ADSL summary ---------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 2
      v Structure (inferred): one record per USUBJID
      Treatment (ARM): Drug A: 1 | Placebo: 1
      Populations: SAFFL 2
      Subject flow: treated 2
      x 1 subject has different `ARM` and `ACTARM`
      v Checks passed: no missing TRTSDT in the safety population; population flags
        only contain Y/N/NA

# print.summary_admiral_df Test 21: BDS formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 5
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISITN
      Analysis visits (AVISIT): BASELINE, WEEK 2, WEEK 4
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          4         2       3        1   60      65   70
        SYSBP          1         1       1        0  120     120  120
      Derived records (DTYPE): LOCF: 1
      v Checks passed: AVISIT and AVISITN are consistent; one PARAM/AVALU per PARAMCD

# print.summary_admiral_df Test 25: OCCDS formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: OCCDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 4
      v Structure (inferred): one record per USUBJID, ASEQ
      Distinct terms: AEBODSYS 2 | AEDECOD 3
      Treatment emergent (TRTEMFL): 3 of 4 records
      Severity/grade (AESEV): MILD: 2 | MODERATE: 1 | SEVERE: 1
      Serious (AESER): 1
      v Checks passed: occurrence flags are unique per subject and level; no missing
        ASTDT; no treatment-emergent record before TRTSDT

