# derive_param_tte Test 8: errors if PARAMCD and by_vars are not one to one

    Code
      derive_param_tte(dataset_adsl = adsl, by_vars = exprs(AEDECOD), start_date = TRTSDT,
      event_conditions = list(ttae), censor_conditions = list(eos), source_datasets = list(
        adsl = adsl, ae = ae), set_values_to = exprs(PARAMCD = "TTAE", PARCAT2 = AEDECOD))
    Condition
      Error in `derive_param_tte()`:
      ! For some values of "PARAMCD" there is more than one value of "AEDECOD"
      i Call `get_one_to_many_dataset()` to get all one-to-many values.

