# error if source dataset is not available

    Code
      derive_var_merged_ef_msrc(adsl, flag_events = list(flag_event(dataset_name = "cm",
        condition = CMCAT == "ANTI-CANCER"), flag_event(dataset_name = "pr")),
      source_datasets = list(cm = cm, pro = pr), by_vars = exprs(USUBJID), new_var = CANCTRFL)
    Condition
      Error in `derive_var_merged_ef_msrc()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: cm and pro
      i  But, `flag_events[[2]]$dataset_name = pr`

