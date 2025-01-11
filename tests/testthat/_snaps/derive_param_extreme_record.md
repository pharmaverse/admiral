# derive_param_extreme_record Test 1: deprecation message if function is called

    Code
      df <- derive_param_extreme_record(dataset = aevent, sources = list(
        records_source(dataset_name = "cm", filter = CMDECOD == "ACT", new_vars = exprs(
          ADT = convert_dtc_to_dt(CMSTDTC), AVALC = CMDECOD)), records_source(
          dataset_name = "pr", filter = PRDECOD == "ACS", new_vars = exprs(ADT = convert_dtc_to_dt(
            PRSTDTC), AVALC = PRDECOD))), source_datasets = list(cm = cm, pr = pr),
      by_vars = exprs(STUDYID, USUBJID), order = exprs(ADT), mode = "first",
      set_values_to = exprs(PARAMCD = "FIRSTACT", PARAM = "First Anti-Cancer Therapy"))
    Message
      `derive_param_extreme_record()` was deprecated in admiral 1.2.0.
      i Please use `derive_extreme_event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation

