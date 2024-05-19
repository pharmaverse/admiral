# derive_var_dthcaus Test 11: error if source dataset is not available

    Code
      derive_var_dthcaus(adsl, source_datasets = list(ae = ae, dd = ds), src_ae,
      src_ds)
    Condition
      Error in `derive_var_dthcaus()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ae and dd
      i  But, `sources[[2]]$dataset_name = ds`

