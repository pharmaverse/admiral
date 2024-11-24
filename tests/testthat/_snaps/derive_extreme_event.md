# derive_extreme_event Test 8: error if source dataset not available

    Code
      derive_extreme_event(adhy, by_vars = exprs(USUBJID), events = list(event(
        dataset_name = "adyh", condition = is.na(CRIT1FL), set_values_to = exprs(
          AVALC = "N")), event(condition = CRIT1FL == "Y", mode = "last",
      set_values_to = exprs(AVALC = "Y"))), source_datasets = list(adhy = adhy),
      tmp_event_nr_var = event_nr, order = exprs(event_nr, AVISITN), mode = "first",
      keep_source_vars = exprs(AVISITN), set_values_to = exprs(PARAMCD = "ALK2",
        PARAM = "ALKPH <= 2 times ULN"))
    Condition
      Error in `derive_extreme_event()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: adhy
      i  But, `events[[1]]$dataset_name = adyh`

