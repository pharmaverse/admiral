# derive_vars_joined Test 8: error if new_vars are already in dataset

    Code
      derive_vars_joined(myd, dataset_add = myd, order = exprs(day), join_type = "all",
      mode = "last", filter_join = day < day.join)
    Condition
      Error in `derive_vars_joined()`:
      ! The variables `day` and `val` in `dataset_add` have naming conflicts with `dataset`, please make the appropriate modifications to `new_vars`.

