# derive_vars_transposed Test 3: filtering the merge dataset works
          with relationship 'many-to-one'

    Code
      derive_vars_transposed(dataset, dataset_merge, by_vars = exprs(USUBJID),
      key_var = TESTCD, value_var = VALUE, filter = TESTCD == "T01", relationship = "many-to-one")
    Output
      # A tibble: 3 x 3
        USUBJID  VAR1   T01
        <chr>   <dbl> <dbl>
      1 P01         3    31
      2 P02        31     3
      3 P03        42    NA

