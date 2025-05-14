# derive_vars_joined_summary Test 3: error if new variable in input dataset

    Code
      derive_vars_joined_summary(dataset = adae, dataset_add = adex, by_vars = exprs(
        USUBJID), filter_join = ADY.join <= ADY, join_type = "all", join_vars = exprs(
        ADY), new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE)))
    Condition
      Error in `derive_vars_joined_summary()`:
      ! The variable `CUMDOSE` in `new_vars` is already in `dataset`
      Please make appropriate modifications to `dataset` or `new_vars`.

