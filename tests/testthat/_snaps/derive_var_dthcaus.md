# derive_var_dthcaus Test 1: deprecation message if function is called

    Code
      src_ae <- dthcaus_source(dataset_name = "ae", filter = AEOUT == "FATAL", date = AEDTHDT,
      mode = "first", dthcaus = AEDECOD)
    Message
      `dthcaus_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      src_ds <- dthcaus_source(dataset_name = "ds", filter = DSDECOD == "DEATH" &
        grepl("DEATH DUE TO", DSTERM), date = convert_dtc_to_dt(DSSTDTC), mode = "first",
      dthcaus = str_to_upper(DSTERM))
    Message
      `dthcaus_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      derive_var_dthcaus(adsl, source_datasets = list(ae = ae, ds = ds), src_ae,
      src_ds)
    Message
      `derive_var_dthcaus()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Output
      # A tibble: 3 x 3
        STUDYID USUBJID DTHCAUS                            
        <chr>   <chr>   <chr>                              
      1 TEST01  PAT01   DEATH DUE TO PROGRESSION OF DISEASE
      2 TEST01  PAT02   <NA>                               
      3 TEST01  PAT03   SUDDEN DEATH                       

# derive_var_dthcaus Test 12: error if source dataset is not available

    Code
      derive_var_dthcaus(adsl, source_datasets = list(ae = ae, dd = ds), src_ae,
      src_ds)
    Condition
      Error in `derive_var_dthcaus()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ae and dd
      i But, `sources[[2]]$dataset_name = ds`

