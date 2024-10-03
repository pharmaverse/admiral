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

