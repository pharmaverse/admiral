# derive_vars_merged Test 5: error if variable in both datasets

    Code
      derive_vars_merged(advs, dataset_add = adsl, by_vars = exprs(USUBJID))
    Condition
      Error in `derive_vars_merged()`:
      ! The variable `STUDYID` is contained in both datasets.
      i Please add them to `by_vars` or remove or rename them in one of the datasets.

# derive_vars_merged Test 11: error if not unique w.r.t the by variables and the order

    Code
      actual <- derive_vars_merged(advs, dataset_add = adsl2, by_vars = exprs(STUDYID,
        USUBJID = ID), order = exprs(ID), mode = "last", check_type = "error",
      duplicate_msg = "Duplicate records present!")
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset contains duplicate records with respect to `STUDYID`, `ID`, and `ID`
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

# derive_vars_merged Test 12: error if variables in missing_values but not in new_vars

    Code
      derive_vars_merged(adsl, dataset_add = advs, by_vars = exprs(USUBJID), order = exprs(
        AVISIT), new_vars = exprs(LASTVIS = str_to_upper(AVISIT)), mode = "last",
      missing_values = exprs(LASTVIS = "UNKNOWN", LASTVISN = -1))
    Condition
      Error in `derive_vars_merged()`:
      ! The variable `LASTVISN` was specified for `missing_values` but not for `new_vars`.

# derive_vars_merged_lookup Test 16: merge lookup table

    Code
      actual <- derive_vars_merged_lookup(vs, dataset_add = param_lookup, by_vars = exprs(
        VSTESTCD, VSTEST), new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE)
    Message
      List of `VSTESTCD` and `VSTEST` not mapped:
      # A tibble: 1 x 2
      VSTESTCD VSTEST
      <chr> <chr>
      1 DIABP Diastolic Blood Pressure
      i Run `admiral::get_not_mapped()` to access the full list.

# derive_vars_merged_lookup Test 18: by_vars with rename

    Code
      actual <- derive_vars_merged_lookup(vs, dataset_add = param_lookup, by_vars = exprs(
        VSTESTCD = TESTCD, VSTEST), new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE)
    Message
      List of `VSTESTCD` and `VSTEST` not mapped:
      # A tibble: 1 x 2
      VSTESTCD VSTEST
      <chr> <chr>
      1 DIABP Diastolic Blood Pressure
      i Run `admiral::get_not_mapped()` to access the full list.

# get_not_mapped Test 19: not all by_vars have records in the lookup table

    Code
      act_vs_param <- derive_vars_merged_lookup(vs, dataset_add = param_lookup,
        by_vars = exprs(VSTESTCD, VSTEST), new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
        print_not_mapped = TRUE)
    Message
      List of `VSTESTCD` and `VSTEST` not mapped:
      # A tibble: 1 x 2
      VSTESTCD VSTEST
      <chr> <chr>
      1 DIABP Diastolic Blood Pressure
      i Run `admiral::get_not_mapped()` to access the full list.

# derive_var_merged_summary Test 24: error if not unique, no order, check_type = error

    Code
      actual <- derive_vars_merged(advs, dataset_add = adsl2, by_vars = exprs(STUDYID,
        USUBJID = ID), order = NULL, check_type = "error")
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset `dataset_add` contains duplicate records with respect to `STUDYID` and `ID`.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

# derive_var_merged_summary Test 25: error if not unique, no order, check_type = warning

    Code
      actual <- derive_vars_merged(advs, dataset_add = adsl2, by_vars = exprs(STUDYID,
        USUBJID = ID), order = NULL, check_type = "warning")
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset `dataset_add` contains duplicate records with respect to `STUDYID` and `ID`.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

# derive_var_merged_summary Test 26: error if not unique, no order, check_type = NULL

    Code
      actual <- derive_vars_merged(advs, dataset_add = adsl2, by_vars = exprs(STUDYID,
        USUBJID = ID), order = NULL, check_type = NULL)
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset `dataset_add` contains duplicate records with respect to `STUDYID` and `ID`.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

