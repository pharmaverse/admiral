# derive_var_extreme_dtm Test 6: error if source dataset is not available

    Code
      derive_var_extreme_dtm(adsl, new_var = LSTALVDTM, source_datasets = list(ea = ae),
      ae_start, mode = "last")
    Condition
      Error in `derive_var_extreme_dtm()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ea
      i  But, `sources[[1]]$dataset_name = ae`

