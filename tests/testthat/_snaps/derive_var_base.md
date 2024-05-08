# derive_var_base Test 4: error if multiple baseline records

    Code
      derive_var_base(input, by_vars = exprs(USUBJID, PARAMCD, BASETYPE), source_var = AVALC,
      new_var = BASEC)
    Condition
      Error in `signal_duplicate_records()`:
      ! Input dataset contains multiple baseline records with respect to `USUBJID`, `PARAMCD`, and `BASETYPE`
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

