# derive_vars_joined Test 8: error if new_vars are already in dataset

    Code
      derive_vars_joined(myd, dataset_add = myd, order = exprs(day), join_type = "all",
      mode = "last", filter_join = day < day.join)
    Condition
      Error in `derive_vars_joined()`:
      ! The variables `day` and `val` in `dataset_add` have naming conflicts with `dataset`, please make the appropriate modifications to `new_vars`.

# derive_vars_joined Test 16: warning if `order` is not unique

    Code
      actual <- derive_vars_joined(adbds, dataset_add = adbds, by_vars = exprs(subj),
      order = exprs(day), new_vars = exprs(prevposval = val), join_type = "before",
      mode = "last", filter_add = val >= 0)
    Condition
      Warning:
      Dataset `dataset` contains duplicate records with respect to `subj` and `day`
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records
      Warning:
      Dataset `dataset_add` contains duplicate records with respect to `subj` and `day`
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records
      Warning:
      Dataset contains duplicate records with respect to `subj`, `tmp_obs_nr_1`, and `day.join`
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

