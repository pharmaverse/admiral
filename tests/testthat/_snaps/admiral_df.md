# set_admiral_attrs Test 9: empty keys are stored, keys not in the dataset are not

    Code
      typo <- set_admiral_attrs(input, keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
    Condition
      Warning:
      Key variable `PARMACD` is not in `dataset`, so no keys were set.

---

    Code
      rekeyed <- set_admiral_attrs(keyed, keys = exprs(USUBJID, PARMACD, AVISIT))
    Condition
      Warning:
      Key variables `PARMACD` and `AVISIT` are not in `dataset`, so no keys were set.
      i The keys set previously are kept: `USUBJID` and `PARAMCD`.

# set_admiral_attrs Test 11: an error is issued for invalid arguments

    Code
      set_admiral_attrs(input, keys = c("USUBJID", "PARAMCD"))
    Condition
      Error in `set_admiral_attrs()`:
      ! Argument `keys` must be a list of variables created by `exprs()`, but is a character vector.
      i Use `chr2vars()` to convert a character vector, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

---

    Code
      set_admiral_attrs(input, dataset_type = "ADAE")
    Condition
      Error in `set_admiral_attrs()`:
      ! Argument `dataset_type` must be equal to one of "ADSL", "BDS", "OCCDS", "TTE", "ADDL", "MDBDS", or "MDOCCDS".

# set_admiral_attrs Test 12: device types are accepted, and checked for SPDEVID

    Code
      actual <- set_admiral_attrs(no_device, dataset_type = "MDBDS")
    Condition
      Warning:
      Dataset type "MDBDS" is a medical device type, but `dataset` does not contain `SPDEVID`.
      i Check that `dataset_type` is correct.

# print.admiral_df Test 25: the name, type and keys are printed in a header

    Code
      print(set_admiral_attrs(input, keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS", dataset_type = "BDS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS, keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_attrs(input, keys = exprs(USUBJID, PARAMCD, AVISIT)))
    Output
      # An admiral dataset (type: BDS (inferred), keys: USUBJID, PARAMCD, AVISIT)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_attrs(input, dataset_type = "BDS"))
    Output
      # An admiral dataset (type: BDS)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_attrs(input, keys = exprs(USUBJID), dataset_type = "TTE"))
    Output
      # An admiral dataset (type: TTE, keys: USUBJID)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_attrs(lookup, keys = exprs(PARAMCD), dataset_name = "PARAMS"))
    Output
      # An admiral dataset (name: PARAMS, keys: PARAMCD)
      # A tibble: 1 x 2
        PARAMCD PARAM                          
        <chr>   <chr>                          
      1 DIABP   Diastolic Blood Pressure (mmHg)

---

    Code
      print(set_admiral_attrs(input, keys = exprs(), dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS (inferred), keys: none)
      # A tibble: 2 x 5
        STUDYID USUBJID PARAMCD AVISIT    AVAL
        <chr>   <chr>   <chr>   <chr>    <dbl>
      1 P       1       DIABP   BASELINE    51
      2 P       1       SYSBP   BASELINE   121

---

    Code
      print(set_admiral_attrs(as.data.frame(input), keys = exprs(USUBJID),
      dataset_name = "ADVS"))
    Output
      # An admiral dataset (name: ADVS, type: BDS (inferred), keys: USUBJID)
        STUDYID USUBJID PARAMCD   AVISIT AVAL
      1       P       1   DIABP BASELINE   51
      2       P       1   SYSBP BASELINE  121

