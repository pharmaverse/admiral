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

# derive_vars_atc Test 5: error if facm not unique

    Code
      derive_vars_atc(dataset = cm, dataset_facm = facm)
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset `dataset_facm` contains duplicate records with respect to `STUDYID`, `USUBJID`, `FAREFID`, and `FATESTCD`
      Please check data and `by_vars` and `id_vars` arguments.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

