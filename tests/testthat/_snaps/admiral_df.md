# summary.admiral_df Test 8: formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 5
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

# print.summary_admiral_df Test 15: the check lines distinguish passed, failed, and not run

    Code
      print(summary(clean))
    Message
      -- clean: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 1
      Observations: 2 | Variables: 6
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISIT
      Analysis visits (AVISIT): BASELINE, WEEK 2
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          2         1       2        0  121   125.5  130
      v Checks passed: at most one baseline per subject and parameter

---

    Code
      print(summary(broken))
    Message
      -- broken: BDS summary ---------------------------------------------------------
      Subjects (USUBJID): 1
      Observations: 2 | Variables: 6
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISIT
      Analysis visits (AVISIT): BASELINE, BASELINE 2
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          2         1       2        0  118   119.5  121
      x 1 subject-parameter combination has more than one baseline record (`ABLFL`)

---

    Code
      print(summary(not_run))
    Message
      -- not_run: BDS summary --------------------------------------------------------
      Subjects (USUBJID): 1
      Observations: 2 | Variables: 5
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISIT
      Analysis visits (AVISIT): BASELINE, WEEK 2
      Parameters (1):
        PARAMCD  records  subjects  visits  missing  min  median  max
        SYSBP          2         1       2        0  121   125.5  130

# print.summary_admiral_df Test 21: BDS formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 5 | Variables: 7
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISITN
      Analysis visits (AVISIT): BASELINE, WEEK 2, WEEK 4
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median  max
        DIABP          4         2       3        1   60      65   70
        SYSBP          1         1       1        0  120     120  120
      Derived records (DTYPE): LOCF: 1

# print.summary_admiral_df Test 30: ADSL comparison formatted output is stable

    Code
      print(summary(input, adsl = adsl))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 3
      Observations: 3 | Variables: 5
      v Structure (inferred): one record per USUBJID, PARAMCD
      Parameters (1):
        PARAMCD  records  subjects  missing  min  median  max
        SYSBP          3         3        0  115     121  130
      Compared with ADSL: 2 of 3 subjects have records (66.7%) | safety population 2
      of 2
      By arm (TRT01A): Active 1/1 | Placebo 1/2
      x 1 subject is not in ADSL (e.g. "9")

# print.summary_admiral_df Test 25: OCCDS formatted output is stable

    Code
      print(summary(input))
    Message
      -- input: OCCDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 4 | Variables: 10
      v Structure (inferred): one record per USUBJID, ASEQ
      Distinct terms: AEBODSYS 2 | AEDECOD 3
      Treatment emergent (TRTEMFL): 3 of 4 records
      Severity/grade (AESEV): MILD: 2 | MODERATE: 1 | SEVERE: 1
      Serious (AESER): 1
      v Checks passed: occurrence flags are unique per subject and level

