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

