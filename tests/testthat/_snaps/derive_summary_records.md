# derive_summary_records Test 6: error if no summary function

    Code
      derive_summary_records(dataset_add = adbds, by_vars = exprs(AVISIT),
      set_values_to = exprs(MEANVIS = AVAL / 2))
    Condition
      Error in `signal_duplicate_records()`:
      ! After summarising, the dataset contains mulitple records with respect to `AVISIT`.
      Please check the `set_values_to` argument if summary functions like `mean()`, `sum()`, ... are used on the right hand side.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

