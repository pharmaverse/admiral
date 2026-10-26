# set_admiral_keys Test 8: empty keys are stored, keys not in the dataset are not

    Code
      typo <- set_admiral_keys(input, keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
    Condition
      Warning:
      Key variable `PARMACD` is not in `dataset`, so no keys were set.

---

    Code
      rekeyed <- set_admiral_keys(keyed, keys = exprs(USUBJID, PARMACD, AVISIT))
    Condition
      Warning:
      Key variables `PARMACD` and `AVISIT` are not in `dataset`, so no keys were set.
      i The keys set previously are kept: `USUBJID` and `PARAMCD`.

# set_admiral_keys Test 10: an error is issued for invalid arguments

    Code
      set_admiral_keys(input, keys = c("USUBJID", "PARAMCD"))
    Condition
      Error in `set_admiral_keys()`:
      ! Argument `keys` must be a list of variables created by `exprs()`, but is a character vector.
      i Use `chr2vars()` to convert a character vector, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

# print.admiral_df Test 34: the name, type and keys are printed in a header

    Code
      print(set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS, keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD, AVISIT)))
    Output
      # An admiral dataset (type: BDS, keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_keys(input, keys = exprs(), dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS, keys: none)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_keys(as.data.frame(input), keys = exprs(USUBJID),
      dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS, keys: USUBJID)
        STUDYID USUBJID PARAMCD   AVISIT AVAL
      1       P       1   DIABP BASELINE   51
      2       P       1   SYSBP BASELINE  121

---

    Code
      print(set_admiral_keys(lookup, keys = exprs(PARAMCD), dataset_name = "PARAMS"))
    Output
      # An admiral dataset (name: PARAMS, keys: PARAMCD)
      # A tibble: 1 x 2
        PARAMCD PARAM                          
        <chr>   <chr>                          
      1 DIABP   Diastolic Blood Pressure (mmHg)

# summary.admiral_df Test 38: duplicate, missing and empty declared keys

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
      ! Structure: not checked, declared key AVISIT is not in the dataset
      Parameters (1):
        PARAMCD  records  subjects  missing  min  median  max
        SYSBP          3         2        0  118     121  130

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

# summary.admiral_df Test 39: the ADSL summary

    Code
      print(actual)
    Message
      -- ADSL summary ----------------------------------------------------------------
      Subjects (USUBJID): 3
      Observations: 3 | Variables: 5
      v Structure: one record per USUBJID
      Subjects per treatment arm (TRT01P): Drug A: 2 | Placebo: 1
      Populations: SAFFL 2 | ITTFL 3

# summary.admiral_df Test 40: the BDS summary

    Code
      print(actual)
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 4 | Variables: 6
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISITN
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
      v Structure (inferred): one record per USUBJID, PARAMCD
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

# summary.admiral_df Test 41: the OCCDS summary

    Code
      print(actual)
    Message
      -- input: OCCDS summary --------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 6
      v Structure (inferred): one record per USUBJID, AESEQ
      Distinct terms: AEBODSYS 2 | AEDECOD 2
      Treatment emergent (TRTEMFL): 2 of 3 records

# summary.admiral_df Test 42: the TTE summary

    Code
      print(actual)
    Message
      -- input: TTE summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 3 | Variables: 6
      v Structure (inferred): one record per USUBJID, PARAMCD
      Events and censored per parameter (CNSR):
        PARAMCD  events  censored
        OS            1         1
        PFS           0         1

# summary.admiral_df Test 43: datasets with little to report do not error

    Code
      print(summary(lookup))
    Message
      -- lookup: Dataset summary -----------------------------------------------------
      Observations: 2 | Variables: 1
      i Structure: could not be determined

---

    Code
      print(summary(empty))
    Message
      -- empty: Dataset summary ------------------------------------------------------
      Subjects (USUBJID): 0
      Observations: 0 | Variables: 2
      i Structure: could not be determined

# print_summary_table Test 45: tables are aligned, counts are not rounded, NA is shown as `.`

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

# summary.admiral_df Test 46: the ADSL subject flow, outcome and demographics

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

# summary.admiral_df Test 47: the BDS per-parameter statistics and derived records

    Code
      print(summary(input))
    Message
      -- input: BDS summary ----------------------------------------------------------
      Subjects (USUBJID): 2
      Observations: 5 | Variables: 7
      v Structure (inferred): one record per USUBJID, PARAMCD, AVISIT, DTYPE
      Parameters (2):
        PARAMCD  records  subjects  visits  missing  min  median    max
        DIABP          1         1       0        1    .       .      .
        SYSBP          4         2       2        1  120     125  130.5
      Analysis visits (AVISIT): BASELINE, WEEK 2
      Derived records (DTYPE): LOCF: 1

