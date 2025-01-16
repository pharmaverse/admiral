# derive_vars_atc Test 2: error if facm not unique

    Code
      derive_vars_atc(dataset = cm, dataset_facm = facm)
    Condition
      Error in `signal_duplicate_records()`:
      ! Dataset `dataset_facm` contains duplicate records with respect to `STUDYID`, `USUBJID`, `FAREFID`, and `FATESTCD`
      Please check data and `by_vars` and `id_vars` arguments.
      i Run `admiral::get_duplicates_dataset()` to access the duplicate records

