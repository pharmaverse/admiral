# set_admiral_keys Test 10: an error is issued for invalid arguments

    Code
      set_admiral_keys(input, keys = c("USUBJID", "PARAMCD"))
    Condition
      Error in `set_admiral_keys()`:
      ! Argument `keys` must be a list of variables created by `exprs()`, but is a character vector.
      i Use `chr2vars()` to convert a character vector, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

# print.admiral_df Test 33: the name, type and keys are printed in a header

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

