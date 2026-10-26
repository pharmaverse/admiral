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

# summary.admiral_df Test 39: the ADSL summary

    Code
      print(actual)
    Output
      # Summary of an admiral dataset: ADSL (ADSL)
      Subjects: 3 | Observations: 3 | Variables: 5
      Record structure: one record per USUBJID (declared keys)
      Subjects per treatment arm (TRT01P): Drug A: 2 | Placebo: 1
      Population flags ("Y"): SAFFL 2 | ITTFL 3

# summary.admiral_df Test 40: the BDS summary

    Code
      print(actual)
    Output
      # Summary of an admiral dataset: input (BDS)
      Subjects: 2 | Observations: 4 | Variables: 6
      Record structure: one record per USUBJID, PARAMCD, AVISITN (inferred keys)
      Parameters (2):
         PARAMCD records subjects
           DIABP       1        1
           SYSBP       3        2
      Analysis visits (2): BASELINE, WEEK 2

# summary.admiral_df Test 41: the OCCDS summary

    Code
      print(actual)
    Output
      # Summary of an admiral dataset: input (OCCDS)
      Subjects: 2 | Observations: 3 | Variables: 6
      Record structure: one record per USUBJID, AESEQ (inferred keys)
      Distinct terms: AEBODSYS 2 | AEDECOD 2
      Treatment-emergent records (TRTEMFL): 2 of 3

# summary.admiral_df Test 42: the TTE summary

    Code
      print(actual)
    Output
      # Summary of an admiral dataset: input (TTE)
      Subjects: 2 | Observations: 3 | Variables: 6
      Record structure: one record per USUBJID, PARAMCD (inferred keys)
      Events and censored per parameter:
         PARAMCD events censored
              OS      1        1
             PFS      0        1

# summary.admiral_df Test 43: datasets with little to report do not error

    Code
      print(summary(lookup))
    Output
      # Summary of an admiral dataset: lookup (other)
      Observations: 2 | Variables: 1
      Record structure: could not be determined

---

    Code
      print(summary(empty))
    Output
      # Summary of an admiral dataset: empty (unknown type)
      Subjects: 0 | Observations: 0 | Variables: 2
      Record structure: could not be determined

