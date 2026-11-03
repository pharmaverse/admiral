# admiral_df Test 9: empty keys are stored, keys not in the dataset are not

    Code
      typo <- admiral_df(input, keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
    Condition
      Warning:
      Key variable `PARMACD` is not in `dataset`, so no keys were set.

---

    Code
      rekeyed <- admiral_df(keyed, keys = exprs(USUBJID, PARMACD, AVISIT))
    Condition
      Warning:
      Key variables `PARMACD` and `AVISIT` are not in `dataset`, so no keys were set.
      i The keys set previously are kept: `USUBJID` and `PARAMCD`.

# admiral_df Test 11: an error is issued for invalid arguments

    Code
      admiral_df(input, keys = c("USUBJID", "PARAMCD"))
    Condition
      Error in `admiral_df()`:
      ! Argument `keys` must be a list of variables created by `exprs()`, but is a character vector.
      i Use `chr2vars()` to convert a character vector, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

---

    Code
      admiral_df(input, dataset_type = "ADAE")
    Condition
      Error in `admiral_df()`:
      ! Argument `dataset_type` must be equal to one of "ADSL", "BDS", "OCCDS", "TTE", "ADDL", "MDBDS", or "MDOCCDS".

# print.admiral_df Test 27: the name, type and keys are printed in a header

    Code
      print(admiral_df(input, keys = exprs(USUBJID, PARAMCD, AVISIT), dataset_name = "ADVS",
      dataset_type = "BDS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS, keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, keys = exprs(USUBJID, PARAMCD, AVISIT)))
    Output
      # An admiral dataset (type: BDS (inferred), keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, dataset_type = "BDS"))
    Output
      # An admiral dataset (type: BDS)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(input, keys = exprs(USUBJID), dataset_type = "TTE"))
    Output
      # An admiral dataset (type: TTE, keys: USUBJID)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(lookup, keys = exprs(PARAMCD), dataset_name = "PARAMS"))
    Output
      # An admiral dataset (name: PARAMS, keys: PARAMCD)
      # A tibble: 1 x 2
        PARAMCD PARAM                          
        <chr>   <chr>                          
      1 DIABP   Diastolic Blood Pressure (mmHg)

---

    Code
      print(admiral_df(input, keys = exprs(), dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS (inferred), keys: none)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(admiral_df(as.data.frame(input), keys = exprs(USUBJID), dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS (inferred), keys: USUBJID)
        STUDYID USUBJID PARAMCD   AVISIT AVAL
      1       P       1   DIABP BASELINE   51
      2       P       1   SYSBP BASELINE  121

