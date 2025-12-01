# derive_vars_transposed Test 3: filter merge dataset 'many-to-one'

    Code
      derive_vars_transposed(dataset, dataset_merge, by_vars = get_admiral_option(
        "subject_keys"), key_var = TESTCD, value_var = VALUE, filter = TESTCD ==
      "T01", relationship = "many-to-one")
    Output
      # A tibble: 3 x 4
        STUDYID USUBJID  VAR1   T01
        <chr>   <chr>   <dbl> <dbl>
      1 STUDY01 P01         3    31
      2 STUDY01 P02        31     3
      3 STUDY01 P03        42    NA

# derive_vars_transposed Test 4: error if `relationship` is unexpected

    Code
      cm %>% derive_vars_transposed(facm, by_vars = exprs(USUBJID, CMREFID = FAREFID),
      id_vars = exprs(FAGRPID), key_var = FATESTCD, value_var = FASTRESC,
      relationship = "one-to-one")
    Condition
      Error in `tryCatch()`:
      ! Each row in `dataset` must match at most 1 row in the transposed `dataset_merge`.
      Row 2 of `dataset` matches multiple rows in the transposed `dataset_merge`.

---

    Code
      cm %>% derive_vars_transposed(facm, by_vars = exprs(USUBJID, CMREFID = FAREFID),
      id_vars = exprs(FAGRPID), key_var = FATESTCD, value_var = FASTRESC,
      relationship = "many-to-one")
    Condition
      Error in `derive_vars_transposed()`:
      ! Each row in `dataset` must match at most 1 row in the transposed `dataset_merge`.
      Row 2 of `dataset` matches multiple rows in the transposed `dataset_merge`.

